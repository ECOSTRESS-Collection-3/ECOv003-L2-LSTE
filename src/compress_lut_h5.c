#include <assert.h>
#include <stdio.h>
#include "fileio.h"

int main(int argc, char *argv[])
{
    int stat;

    // Compression parameters
    const int chunk_factor = 12;
    const int compression_deflate = 9;
    const int n_lines = 576;
    const int n_pixels = 361;
    // Specify compression for the output product
    set_hdf5_2d_compression_on(
        n_lines / chunk_factor, n_pixels / chunk_factor, compression_deflate);

    unsigned int nn;
    unsigned int hh;
    unsigned int mm;
    
    char tmpBff[PATH_MAX];

    hid_t infid = open_hdf5(argv[1], H5F_ACC_RDONLY);
    hid_t outfid =  create_hdf5(argv[2]);
    hid_t data_group = H5Gcreate2(outfid, "/Data", 
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(data_group >= 0);
    hid_t geo_group = H5Gcreate2(outfid, "/Geolocation", 
        H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    assert(geo_group >= 0);
    
    size_t check_4_truncated = -1;

    for (nn = 1; nn <=3; nn++)
    {
        for (hh = 0; hh < 24;  hh += 6)
        {
            for (mm = 1; mm <= 12; mm++)
            {
                char dset_name[256];
                snprintf(dset_name, sizeof(dset_name), "/Data/LUT_cloudBT%d_%02d_%02d", nn, hh, mm);
                Matrix lut_data;
                mat_init(&lut_data);
                stat = hdf5read_mat(infid, dset_name, &lut_data);
                if (stat < 0)
                {                    
                    memset(tmpBff, 0, sizeof(tmpBff));
                    check_4_truncated = snprintf(tmpBff, sizeof(tmpBff), "READ error on dataset %s\n", dset_name);
                    if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(tmpBff))
                    {
                        ERREXIT(113, "Buffer is too small.", NONE);
                    }
                    printf(tmpBff);
                }
                assert(stat >= 0);
                stat = hdf5write_mat(outfid, dset_name, &lut_data);
                if (stat < 0)
                {                    
                    memset(tmpBff, 0, sizeof(tmpBff));
                    check_4_truncated = snprintf(tmpBff, sizeof(tmpBff), "WRITE error on dataset %s\n", dset_name);
                    if (check_4_truncated < 0 || (unsigned) check_4_truncated >= sizeof(tmpBff))
                    {
                        ERREXIT(113, "Buffer is too small.", NONE);
                    }
                    printf(tmpBff);
                }
                assert(stat >= 0);
                mat_clear(&lut_data);
            }
        }
    }
    Matrix lat;
    mat_init(&lat);
    stat = hdf5read_mat(infid, "/Geolocation/Latitude", &lat);
    assert(stat >= 0);
    stat = hdf5write_mat(outfid, "/Geolocation/Latitude", &lat);
    assert(stat >= 0);
    mat_clear(&lat);
    Matrix lon;
    mat_init(&lon);
    stat = hdf5read_mat(infid, "/Geolocation/Longitude", &lon);
    assert(stat >= 0);
    stat = hdf5write_mat(outfid, "/Geolocation/Longitude", &lon);
    assert(stat >= 0);
    mat_clear(&lon);
    close_hdf5(infid);
    close_hdf5(outfid);
    return 0;
}
