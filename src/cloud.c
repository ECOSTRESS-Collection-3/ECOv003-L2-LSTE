// See details in cloud.h
//
#include "lste_lib.h"
#include "cloud.h"
#include "error.h"
#include "interps.h"
#include "matrix.h"
#include "metadata.h"
#include "smooth2d.h"
#include "tes_util.h"
#include "fileio.h"

#include <string.h>
#include <stdio.h>

// #define DEBUG

#ifdef DEBUG
    int info;
    hid_t dbg_id;
    char debug_output[PATH_MAX];
    char debugf[PATH_MAX];
    char debugfn[PATH_MAX];    
    char tmpBff[PATH_MAX];
#endif

void process_cloud(RAD *rad,
                   const char* bt11_lut_file,
                   double* rad_lut_data[],
                   int rad_lut_nlines,
                   const char* cloud_filename,
                   const char* product_path,
                   int collection,
                   Mat3d *emis)
{
    int check_4_truncated = -1;
    #ifdef DEBUG                
        check_4_truncated = snprintf(debug_output, sizeof(debug_output), "%s", product_path);
        if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(debug_output))
        {
            ERREXIT(113, "Buffer is too small.", NONE);
        }
    #endif
    const unsigned int BAND_4 = rad->nchannels - 2;    
    const unsigned int COL_1 = 0;
    const unsigned int COL_5 = 4;   
    const double slapse = 6.5;

    // Convert radiances to brightness temperature
    unsigned int nrows = rad->Rad[0].size1;
    unsigned int ncols = rad->Rad[0].size2;
    Matrix TB4; mat_init(&TB4);
    mat_alloc(&TB4, nrows, ncols);    
    unsigned int i;
    double v;
    for (i = 0; i < nrows * ncols; i++)
    {
        v = interp1d_npts(rad_lut_data[COL_5], rad_lut_data[COL_1],
            rad->Rad[BAND_4].vals[i], rad_lut_nlines);
        TB4.vals[i] = v;
    }

    // If collection == 3, compute Emismm.
    Matrix Emismm; mat_init(&Emismm);
    if (collection == 3)
    {
        mat_alloc(&Emismm, nrows, ncols);
        int r, c;
        for (r = 0; r < nrows; r++)
        {
            for (c = 0; c < ncols; c++)
            {
                mat_set(&Emismm, r, c, 
                       (mat3d_get(emis, 3, r, c) + mat3d_get(emis, 4, r, c)) / 2.0);
#ifdef DEBUG_ROW
                if (r == DEBUG_ROW && c == DEBUG_COL)
                {
                    printf("Emis4[DEBUG_ROW,DEBUG_COL]=%.4f\n", mat3d_get(emis, 3, r, c));
                    printf("Emis5[DEBUG_ROW,DEBUG_COL]=%.4f\n", mat3d_get(emis, 4, r, c));
                    printf("Emism[DEBUG_ROW,DEBUG_COL]=%.4f\n", mat_get(&Emismm, r, c));
                }
#endif
            }
        }
        // Smooth with 5x5 moving window. Window size is (2*Nr + 1)x(2*Nc+1).
        smooth2d(Emismm.vals, nrows, ncols, 2, 2);
#ifdef DEBUG_ROW
        printf("Emismm[DEBUG_ROW,DEBUG_COL]=%.4f\n", mat_get(&Emismm, 3719, 4291));
#endif
    }

    // Cloud test 1: Brightness temperature look-up-table (LUT) approach
    // There are now 4 LUT files, one for each 6 hours.
    char* hour_string[] = {"00", "06", "12", "18"};
    float monthdays[] = {0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    hid_t lutid[4];
    char lut_name[256];    
    check_4_truncated = snprintf(lut_name, sizeof(lut_name), "%s", bt11_lut_file);
    if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(lut_name))
    {
        ERREXIT(113, "Buffer is too small.", NONE);
    } 
    char *hour_loc = strstr(lut_name, "??");
    if (hour_loc == NULL)
        ERREXIT(14, "cloudLUT PGE Parameter is missing \"??\" for hour: %s", bt11_lut_file);
    for (i = 0; i < 4; i++)
    {
        hour_loc[0] = hour_string[i][0];
        hour_loc[1] = hour_string[i][1];
        lutid[i] = open_hdf5(lut_name, H5F_ACC_RDONLY);
        #ifdef DEBUG            
            memset(tmpBff, 0, sizeof(tmpBff));
            check_4_truncated = snprintf(tmpBff, sizeof(tmpBff), "lutid[%d] = %s\n", i, lut_name);
            if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(tmpBff))
            {
                ERREXIT(113, "Buffer is too small.", NONE);
            }
            printf(tmpBff);
        #endif //DEBUG
        if (lutid[i] < 0)
            ERREXIT(14, "Unable to open LUT: %s", lut_name);
    }
    int stat = -1;

    Matrix t_lut_lat; mat_init(&t_lut_lat);
    stat = hdf5read_mat(lutid[3], "/Geolocation/Latitude", &t_lut_lat);
    if (stat < 0)
        ERREXIT(14, "Error reading /Geolocation/Latitude from %s", lut_name);
    Matrix lut_lat; mat_init(&lut_lat);
    mat_transpose(&lut_lat, &t_lut_lat);
    mat_clear(&t_lut_lat);

    Matrix t_lut_lon; mat_init(&t_lut_lon);
    stat = hdf5read_mat(lutid[3], "/Geolocation/Longitude", &t_lut_lon);
    if (stat < 0)
        ERREXIT(14, "Error reading /Geolocation/Longitude from %s", lut_name);
    Matrix lut_lon; mat_init(&lut_lon);
    mat_transpose(&lut_lon, &t_lut_lon);
    mat_clear(&t_lut_lon);

    // Get ECOSTRESS observation time info    
    const char* start_pos = rindex(cloud_filename, '/');
    if (start_pos == NULL)
        start_pos = cloud_filename;
    const char* L2_CLOUD = strstr(start_pos, "_CLOUD");
    if (L2_CLOUD == NULL)
        ERREXIT(14, "_CLOUD not in filename %s", cloud_filename);
    int hh, mm, mth, day;
    sscanf(L2_CLOUD+21, "%2d%2d", &mth, &day);
    sscanf(L2_CLOUD+26, "%2d%2d", &hh, &mm);

    float hrfrac = hh + (float)mm/60.0;

    // Find nearest Cloud LUT UTC time
    unsigned int Ztime0 = 6 * (hh / 6);
    unsigned int Ztime1 = (Ztime0 + 6) % 24;
    char cloudvar1[32];
    char cloudvar2[32];

    Matrix BTgrid1; mat_init(&BTgrid1);
    mat_alloc(&BTgrid1, nrows,  ncols);
    Matrix BTgrid2; mat_init(&BTgrid2);
    mat_alloc(&BTgrid2, nrows,  ncols);
    Matrix BTout[4];
    for (i = 1; i <= 3; i++)
    {
        mat_init(&BTout[i]);
        mat_alloc(&BTout[i], nrows, ncols);
    }

    // Find bracketing months to current ECOSTRESS observation month
    float obsmthfrac = (float)mth + (float)day/monthdays[mth] - 0.5; // observation month fraction. LUT month represents middle of the month
    int mth1 = (int)obsmthfrac;
    int mth1_1_12 = mth1;
    if (mth1 < 1) mth1_1_12 = 12;
    int mth2 = mth1 + 1;
    if (mth2 > 12) mth2 = 1;
    Matrix BT1; mat_init(&BT1);
    Matrix t_BT1; mat_init(&t_BT1);
    Matrix BT2; mat_init(&BT2);
    Matrix t_BT2; mat_init(&t_BT2);
    unsigned int lut0 = Ztime0 / 6;
    unsigned int lut1 = Ztime1 / 6;

    int LUTthresh;
    for (LUTthresh = 1; LUTthresh <= 3; LUTthresh++)
    {        
        // First month --------------------------------------------------------
        check_4_truncated = snprintf(cloudvar1, sizeof(cloudvar1), "/Data/LUT_cloudBT%d_%02d_%02d", 
                LUTthresh, Ztime0, mth1_1_12);
        if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(cloudvar1))
        {
            ERREXIT(113, "Buffer is too small.", NONE);
        } 
        check_4_truncated = snprintf(cloudvar2, sizeof(cloudvar2), "/Data/LUT_cloudBT%d_%02d_%02d", 
                LUTthresh, Ztime1, mth1_1_12);
        if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(cloudvar2))
        {
            ERREXIT(113, "Buffer is too small.", NONE);
        } 
        #ifdef DEBUG
            printf("First Month: ");            
            memset(tmpBff, 0, sizeof(tmpBff));
            check_4_truncated = snprintf(tmpBff, sizeof(tmpBff), "cloudvar1 = %s\n", cloudvar1);
            if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(tmpBff))
            {
                ERREXIT(113, "Buffer is too small.", NONE);
            }
            printf(tmpBff);           
            check_4_truncated = snprintf(tmpBff, sizeof(tmpBff), "cloudvar2 = %s\n", cloudvar2);
            if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(tmpBff))
            {
                ERREXIT(113, "Buffer is too small.", NONE);
            }
            printf(tmpBff);
        #endif //DEBUG
        stat = hdf5read_mat(lutid[lut0], cloudvar1, &t_BT1);
        if (stat < 0)
        {
            hour_loc[0] = hour_string[lut0][0];
            hour_loc[1] = hour_string[lut0][1];
            ERREXIT(14, "Error reading %s from %s", cloudvar1, lut_name);
        }
        mat_transpose(&BT1, &t_BT1);
        mat_clear(&t_BT1);

        stat = hdf5read_mat(lutid[lut1], cloudvar2, &t_BT2);
        if (stat < 0)
        {
            hour_loc[0] = hour_string[lut1][0];
            hour_loc[1] = hour_string[lut1][1];
            ERREXIT(14, "Error reading %s from %s", cloudvar2, lut_name);
        }
        mat_transpose(&BT2, &t_BT2);               
        mat_clear(&t_BT2);

        // Grid thresholds onto ECOSTRESS scene
        double *inds[2] = {BT1.vals, BT2.vals};
        double *outds[2] = {BTgrid1.vals, BTgrid2.vals};
        multi_interp2(lut_lat.vals, lut_lon.vals, lut_lat.size1, lut_lat.size2,
            rad->Lat.vals, rad->Lon.vals, rad->Lat.size1 * rad->Lat.size2,
            inds, outds, 2);

        // interpolate between bracketing hours
        double m = (hrfrac - Ztime0) / 6.0;
        Matrix BTm1; mat_init(&BTm1);
        mat_alloc(&BTm1, nrows, ncols);
        for (i = 0; i < nrows * ncols; i++)
        {
            BTm1.vals[i] = BTgrid1.vals[i] + m * (BTgrid2.vals[i] - BTgrid1.vals[i]);
        }
                
        // Second month -------------------------------------------------------
        check_4_truncated = snprintf(cloudvar1, sizeof(cloudvar1), "/Data/LUT_cloudBT%d_%02d_%02d", 
                LUTthresh, Ztime0, mth2);
        if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(cloudvar1))
        {
            ERREXIT(113, "Buffer is too small.", NONE);
        } 
        check_4_truncated = snprintf(cloudvar2, sizeof(cloudvar2),"/Data/LUT_cloudBT%d_%02d_%02d", 
                LUTthresh, Ztime1, mth2);
        if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(cloudvar2))
        {
            ERREXIT(113, "Buffer is too small.", NONE);
        } 
        #ifdef DEBUG
            printf("Second Month: ");
            memset(tmpBff, 0, sizeof(tmpBff));            
            check_4_truncated = snprintf(tmpBff, sizeof(tmpBff), "cloudvar1 = %s\n", cloudvar1);
            if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(tmpBff))
            {
                ERREXIT(113, "Buffer is too small.", NONE);
            }
            printf(tmpBff);            
            check_4_truncated = snprintf(tmpBff, sizeof(tmpBff), "cloudvar2 = %s\n", cloudvar2);
            if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(tmpBff))
            {
                ERREXIT(113, "Buffer is too small.", NONE);
            }
            printf(tmpBff);
        #endif //DEBUG
        stat = hdf5read_mat(lutid[lut0], cloudvar1, &t_BT1);
        if (stat < 0)
        {
            hour_loc[0] = hour_string[lut0][0];
            hour_loc[1] = hour_string[lut0][1];
            ERREXIT(14, "Error reading %s from %s", cloudvar1, lut_name);
        }
        mat_transpose(&BT1, &t_BT1);                           
        mat_clear(&t_BT1);

        stat = hdf5read_mat(lutid[lut1], cloudvar2, &t_BT2);
        if (stat < 0)
        {
            hour_loc[0] = hour_string[lut1][0];
            hour_loc[1] = hour_string[lut1][1];
            ERREXIT(14, "Error reading %s from %s", cloudvar2, lut_name);
        }
        mat_transpose(&BT2, &t_BT2);        
        mat_clear(&t_BT2);

        multi_interp2(lut_lat.vals, lut_lon.vals, lut_lat.size1, lut_lat.size2,
            rad->Lat.vals, rad->Lon.vals, rad->Lat.size1 * rad->Lat.size2,
            inds, outds, 2);                     

        // interpolate between bracketing hours
        m = (hrfrac - Ztime0) / 6.0;
        Matrix BTm2; mat_init(&BTm2);
        mat_alloc(&BTm2, nrows, ncols);
        for (i = 0; i < nrows * ncols; i++)
        {
            BTm2.vals[i] = BTgrid1.vals[i] + m * (BTgrid2.vals[i] - BTgrid1.vals[i]);
        }
        // Interpolation between the two months above -------------------------
        double mfobs = obsmthfrac - (double)mth1;
        double BTgridmdiff;
        double BTgridmfrac;
        for (i = 0; i < nrows * ncols; i++)
        {
            BTgridmdiff = BTm2.vals[i] - BTm1.vals[i];
            BTgridmfrac = mfobs * BTgridmdiff;
            BTout[LUTthresh].vals[i] = BTm1.vals[i] + BTgridmfrac;
        }

        // Detected clouds have lower brightness temperatures than thresholds
        // modified by lapse rate.
        double tba;
        for (i = 0; i < nrows * ncols; i++)
        {
            tba = rad->El.vals[i] * slapse;
            if (isnan(TB4.vals[i]))
                BTout[LUTthresh].vals[i] = 255;
            else
                BTout[LUTthresh].vals[i] = BTout[LUTthresh].vals[i] - tba;
        }
        mat_clear(&BT1);
        mat_clear(&BT2);
        
        #ifdef DEBUG            
            memset(debugfn, 0, sizeof(debugfn));
            memset(debugf, 0, sizeof(debugf));
        
            check_4_truncated = snprintf(debugfn, sizeof(debugfn), "/BTout_%d.h5", LUTthresh);
            if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(debugfn))
            {
                ERREXIT(113, "Buffer is too small.", NONE);
            }  
            check_4_truncated = snprintf(debugf, sizeof(debugf), "%s/BTout_%d.h5",debug_output,LUTthresh); 
            if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(debugf))
            {
                ERREXIT(113, "Buffer is too small.", NONE);
            }            
            dbg_id = create_hdf5(debugf);            
            if (dbg_id < 0)
            {                
                memset(tmpBff, 0, sizeof(tmpBff));
                check_4_truncated = snprintf(tmpBff, sizeof(tmpBff), "\nWARNING: Could not open %s\n",debugfn);
                if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(tmpBff))
                {
                    ERREXIT(113, "Buffer is too small.", NONE);
                }
                printf(tmpBff);
            }
            
            info = hdf5write_mat(dbg_id, "/BTout", &BTout[LUTthresh]);
            if (info < 0)
            {                
                memset(tmpBff, 0, sizeof(tmpBff));
                check_4_truncated = snprintf(tmpBff, sizeof(tmpBff), "\nWARNING: Could not write BTout%d to %s\n", LUTthresh, debugfn);
                if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(tmpBff))
                {
                    ERREXIT(113, "Buffer is too small.", NONE);
                }
                printf(tmpBff);
            }
            close_hdf5(dbg_id);
            
            memset(debugfn, 0, sizeof(debugfn));
            memset(debugf, 0, sizeof(debugf));            
        #endif
    }

    MatUint8 cloud1;
    mat_uint8_init(&cloud1);
    mat_uint8_alloc(&cloud1, nrows, ncols);
    MatUint8 cloudconf;
    mat_uint8_init(&cloudconf);
    mat_uint8_alloc(&cloudconf, nrows, ncols);
    MatUint8 emis_cloud;
    mat_uint8_init(&emis_cloud);
    mat_uint8_alloc(&emis_cloud, nrows, ncols);
    uint8 confidentcloud, probablycloud, probablyclear, confidentclear;
#ifdef DEBUG_ROW
    int dix = DEBUG_ROW * 8582 + DEBUG_COL;
#endif
    for (i = 0; i < nrows * ncols;  i++)
    {
        confidentcloud = TB4.vals[i] <= BTout[1].vals[i];
        probablycloud  = (BTout[1].vals[i] < TB4.vals[i]) && 
                         (TB4.vals[i] <= BTout[2].vals[i]);
        probablyclear  = (BTout[2].vals[i] < TB4.vals[i]) && 
                         (TB4.vals[i] <= BTout[3].vals[i]);
        confidentclear = TB4.vals[i] > BTout[3].vals[i];

#ifdef DEBUG_ROW
        if (i == dix)
        {
            printf("Emismm %.4f", Emismm.vals[i]);
            if (confidentcloud)printf(" confident cloud\n");
            if (probablycloud)printf(" probably cloud\n");
            if (probablyclear)printf(" probably clear\n");
            if (confidentclear)printf(" confident clear\n");
        }
#endif
        if (collection == 3)
        {
            if (Emismm.vals[i] < 0.92)
            {
                // Indicate emis_cloud results
                emis_cloud.vals[i] = 1;
            }
            else
            {
                emis_cloud.vals[i] = 0;
            }
        }

        if (isnan(TB4.vals[i]))
        {
            cloud1.vals[i] = 255;
            cloudconf.vals[i] = 255;
            emis_cloud.vals[i] = 255;
        }
        else if (confidentclear)
        {
            cloud1.vals[i] = 0;
            cloudconf.vals[i] = 0;
        }
        else if (probablyclear)
        {
            cloud1.vals[i] = 0;
            cloudconf.vals[i] = 1;
        }
        else if (probablycloud)
        {
            
            // At altitude > 2 Km consider probablycloud as clear.
            if (rad->El.vals[i] > 2.0)
                cloud1.vals[i] = 0;
            else
                cloud1.vals[i] = 1;
        }
        else if (confidentcloud)
        {
            cloud1.vals[i] = 1;
            cloudconf.vals[i] = 3;
        }
        else
        {
            cloud1.vals[i] = 255;
            cloudconf.vals[i] = 255;
        }
    }        

    // Cleanup cloud test 1
    mat_clear(&lut_lat);
    mat_clear(&lut_lon);
    mat_clear(&BTgrid1);
    mat_clear(&BTgrid2);
    for (i = 0; i < 4; i++)
    {
        close_hdf5(lutid[i]);
        if (i > 0)
            mat_clear(&BTout[i]);
    }

    int nquery_points = nrows * ncols;
    /* Skipping test 2 now.

    // Cloud test 2

    // Read Cloud BTdiff data
    Vector T11; vec_init(&T11);
    Matrix Tmax; mat_init(&Tmax);
    Vector VA;  vec_init(&VA);
    h5fid = open_hdf5(btdiff_file, H5F_ACC_RDONLY);
    if (h5fid == -1)
    {
        ERREXIT(14, "Unable to open cloud BTdiff file %s", btdiff_file);
    }
    stat = hdf5read_vec(h5fid, "/SDS/T11", &T11);
    if (stat == -1)
    {
        ERREXIT(14, "Unable to read /SDS/T11 from %d", btdiff_file);
    }
    stat = hdf5read_mat(h5fid, "/SDS/Tmax", &Tmax);
    if (stat == -1)
    {
        ERREXIT(14, "Unable to read /SDS/Tmax from %d", btdiff_file);
    }
    stat = hdf5read_vec(h5fid, "/SDS/VA", &VA);
    if (stat == -1)
    {
        ERREXIT(14, "Unable to read /SDS/VA from %d", btdiff_file);
    }
    close_hdf5(h5fid);

    // [1.1.2] Reduce Tmax by 10%
    int iTmax;
    for (iTmax = 0; iTmax < Tmax.size1 * Tmax.size2; iTmax++)
    {
        Tmax.vals[iTmax] *= 0.9;
    }

    // Set up BT diff meshgrid for interpolants
    int btdiff_nrows = T11.size - 1;
    int btdiff_ncols = VA.size;
    Matrix btdiff_X; mat_init(&btdiff_X);
    mat_alloc(&btdiff_X, btdiff_nrows, btdiff_ncols);
    Matrix btdiff_Y; mat_init(&btdiff_Y);
    mat_alloc(&btdiff_Y, btdiff_nrows, btdiff_ncols);
    int btrow, btcol;
    for (btrow = 0; btrow < btdiff_nrows; btrow++)
    {
        for (btcol = 0; btcol < btdiff_ncols; btcol++)
        {
            mat_set(&btdiff_X, btrow, btcol, T11.vals[btrow+1]);
            mat_set(&btdiff_Y, btrow, btcol, VA.vals[btcol]);
        }
    }

    Matrix TB4temp; mat_init(&TB4temp);
    mat_alloc(&TB4temp, nrows, ncols);
    Matrix SZtemp; mat_init(&SZtemp);
    mat_alloc(&SZtemp, nrows, ncols);
    for (i = 0; i < nquery_points; i++)
    {
        v = TB4.vals[i];
        if (v < 260.0) v = 260.0;
        TB4temp.vals[i] = v;
        v = rad->Satzen.vals[i];
        if (v < 3.0)
            v = 3.0;
        SZtemp.vals[i] = v;
    }

    // Interpolate Tmax onto granule based on TB4 and SZ
    Matrix btdifft_mat; mat_init(&btdifft_mat);
    mat_alloc(&btdifft_mat, nrows, ncols);
    interp2(btdiff_X.vals, btdiff_Y.vals, Tmax.vals, btdiff_nrows, btdiff_ncols,
            TB4temp.vals, SZtemp.vals, btdifft_mat.vals, nquery_points);

    double Tdiff1, btdifft1;
    MatUint8 cloud2;
    mat_uint8_init(&cloud2);
    mat_uint8_alloc(&cloud2, nrows, ncols);
    for (i = 0; i < nquery_points; i++)
    {
        Tdiff1 = TB4.vals[i] - TB5.vals[i];
        // Increase threshold above 1km elevation to account for
        // elevation dependent anomalies.
        btdifft1 = btdifft_mat.vals[i] + 1.0;
        if (isnan(Tdiff1) || isnan(btdifft1))
            cloud2.vals[i] = 255;
        else
            cloud2.vals[i] = Tdiff1 > btdifft1;
    }

    // Contiguity  - set cloudy pixels with 3 or less cloudy neighbors to clear
    const int NumNeighbors = 3;
    constrict_cloud(&cloud2, NumNeighbors);
    // Repeat a second time.
    constrict_cloud(&cloud2, NumNeighbors);
    */

    // Final cloud
    MatUint8 cloud_final;
    mat_uint8_init(&cloud_final);
    mat_uint8_alloc(&cloud_final, nrows, ncols);
    for (i = 0; i < nquery_points; i++)
    {        
        cloud_final.vals[i] = cloud1.vals[i];
    }

    // Cleanup intermediate matrices.
    mat_clear(&TB4);    
    mat_clear(&Emismm);

    // Now outputting cloud1 as Cloud_final and cloudconf as Cloud_confidence.
 
    // Output the Cloud_test_1 dataset.
    hid_t cloudout = create_hdf5(cloud_filename);
    if (cloudout < 0)
        ERREXIT(16, "Could not create output file %s", cloud_filename);
    
    hid_t group[4];
    char group_name[PATH_MAX];    
    check_4_truncated = snprintf(group_name, sizeof(group_name), "%s", "/HDFEOS"); 
    if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(group_name))
    {
        ERREXIT(113, "Buffer is too small.", NONE);
    }
    group[0] = H5Gcreate2(cloudout, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group[0] == -1)
        ERREXIT(19, "Could not create group %s in %s", group_name, cloud_filename);   
    memset(group_name, 0, sizeof(group_name));
    check_4_truncated = snprintf(group_name, sizeof(group_name), "%s", "/HDFEOS/GRIDS");
    if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(group_name))
    {
        ERREXIT(113, "Buffer is too small.", NONE);
    }
    group[1] = H5Gcreate2(cloudout, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group[1] == -1)
        ERREXIT(19, "Could not create group %s in %s", group_name, cloud_filename);    
    memset(group_name, 0, sizeof(group_name));
    check_4_truncated = snprintf(group_name, sizeof(group_name), "%s", "/HDFEOS/GRIDS/ECO_L2G_CLOUD_70m");
    if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(group_name))
    {
        ERREXIT(113, "Buffer is too small.", NONE);
    }
    group[2] = H5Gcreate2(cloudout, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group[2] == -1)
        ERREXIT(19, "Could not create group %s in %s", group_name, cloud_filename);    
    memset(group_name, 0, sizeof(group_name));
    check_4_truncated = snprintf(group_name, sizeof(group_name), "%s", "/HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields");
    if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(group_name))
    {
        ERREXIT(113, "Buffer is too small.", NONE);
    }
    group[3] = H5Gcreate2(cloudout, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group[3] == -1)
        ERREXIT(19, "Could not create group %s in %s", group_name, cloud_filename);

    stat = hdf5write_mat_uint8(cloudout, "/HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields/Cloud_confidence", &cloudconf);
    if (stat == -1)
        ERREXIT(17, "Error writing /HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields/Cloud_confidence to %s", cloud_filename);

    // Add attributes to the output dataset.
    uint8 minCM = 0;
    uint8 maxCM = 1;
    uint8 minCONF = 0;
    uint8 maxCONF = 3;
    uint8 fillCM = 255;
    stat = writeDatasetMetadataHdf5(cloudout, "/HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields/Cloud_confidence",
            "Cloud confidence",
            "0=ConfidentClear\n1=ProbablyClear\n2=ProbablyCloud\n3=ConfidentCloud",
            DFNT_UINT8, &minCONF, &maxCONF,
            DFNT_UINT8, &fillCM, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(17, "Failed to write attributes to dataset "
            "/HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields/Cloud_confidence in %s",
            cloud_filename);

    /* Skipping cloud2.
    // Output Cloud_test_2.
    stat = hdf5write_mat_uint8(cloudout, "/SDS/Cloud_test_2", &cloud2);
    if (stat == -1)
        ERREXIT(17, "Error writing /SDS/Cloud_test_2 to %s", cloud_filename);

    // Add attributes to the output dataset.
    stat = writeDatasetMetadataHdf5(cloudout, "/SDS/Cloud_test_2",
            "Brightness temperature difference test",
            "1 = cloudy",
            DFNT_UINT8, &minCM, &maxCM,
            DFNT_UINT8, &fillCM, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(17, "Failed to write attributes to dataset "
            "/SDS/Cloud_test_2 in %s",
            cloud_filename);
    */

    // Output Cloud_test_final. Use data from cloud1.    
    stat = hdf5write_mat_uint8(cloudout, "/HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields/Cloud_final", &cloud1);
    if (stat == -1)
        ERREXIT(17, "Error writing /HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields/Cloud_final to %s", cloud_filename);

    // Add attributes to the output dataset.
    stat = writeDatasetMetadataHdf5(cloudout, "/HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields/Cloud_final",
            "Final cloud mask",
            "1 = cloudy",
            DFNT_UINT8, &minCM, &maxCM,
            DFNT_UINT8, &fillCM, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(17, "Failed to write attributes to dataset "
            "/HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields/Cloud_final in %s",
            cloud_filename);

    // Output emis_cloud. Use data from emis_cloud.
    stat = hdf5write_mat_uint8(cloudout, "/HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields/emis_cloud", &emis_cloud);
    if (stat == -1)
        ERREXIT(17, "Error writing /HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields/emis_cloud to %s", cloud_filename);

    // Add attributes to the output dataset.
    stat = writeDatasetMetadataHdf5(cloudout, "/HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields/emis_cloud",
            "emissivity detected cloud",
            "1 = cloudy",
            DFNT_UINT8, &minCM, &maxCM,
            DFNT_UINT8, &fillCM, 1.0, 0.0);
    if (stat < 0)
        ERREXIT(17, "Failed to write attributes to dataset "
            "/HDFEOS/GRIDS/ECO_L2G_CLOUD_70m/Data Fields/emis_cloud in %s",
            cloud_filename);

    close_hdf5(cloudout);

    // Cleanup data matrices.
    mat_uint8_clear(&cloud1);
    //mat_uint8_clear(&cloud2);
    //mat_uint8_clear(&cloud_final);
    mat_uint8_clear(&cloudconf);
    mat_uint8_clear(&emis_cloud);
}
