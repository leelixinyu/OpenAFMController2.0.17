
#ifdef _MSC_VER
#define GWYSCANDLL_API __declspec(dllexport) 
#else
#define GWYSCANDLL_API
#endif


#include <stdbool.h>

/*!\defgroup GwyscanSave Save Gwyddion File
* @{
*/
/**@}*/

/*!\defgroup GwyscanCreate Create Scan Path
* @{
*/
/**@}*/

/*!\defgroup GwyscanRefine Refine Scan Path
* @{
This section is under construction. 
*/
/**@}*/

/*!\defgroup GwyscanModify Modify Scan Path
* @{
*/
/**@}*/


/*!\defgroup GwyscanInterpolation Preview data in regular matrix
* @{
*/
/**@}*/

/*!\defgroup GwyscanEnumerations: Direction, Error
* @{
*/
typedef enum {
   GWYSCAN_DIRECTION_FORWARD = 0,
   GWYSCAN_DIRECTION_REVERSE = 1,
   GWYSCAN_DIRECTION_BOTH = 2,
   GWYSCAN_DIRECTION_BIDIRECTIONAL = 3,
   GWYSCAN_DIRECTION_FORWARD_REFINE = 4,
   GWYSCAN_DIRECTION_FORWARD_ALTERNATEBLOCK = 5,
} GwyscanDirection;

typedef enum {
   GWYSCAN_SUCCESS = 0,
   GWYSCAN_DATA_ERROR = 1,
   GWYSCAN_WRITE_FAILED = 2,
   GWYSCAN_ALLOCATION_FAILED = 3
} GwyscanResult;
/**@}*/


/*!\ingroup GwyscanSave
* @{
*/
GWYSCANDLL_API GwyscanResult gwyscan_save_gwyddion_array(const double *data, int xres, int yres, double xreal, double yreal, const char *xyunit, const char *zunit, const char *description, const char *filename);

GWYSCANDLL_API GwyscanResult gwyscan_add_gwyddion_array(const double *data, int xres, int yres, double xreal, double yreal, const char *xyunit, const char *zunit, const char *description, const char *filename);

GWYSCANDLL_API GwyscanResult gwyscan_save_gwyddion_arrays_general(const double **data, int nchannels, const int *xres, const int *yres, const double *xreal, const double *yreal, const char **xyunit, const char **zunit, const char **description, const char *filename);

GWYSCANDLL_API GwyscanResult gwyscan_save_gwyddion_arrays(const double *data, int nchannels, int xres, int yres, double xreal, double yreal, const char **xyunit, const char **zunit, const char **description, const char *filename);

GWYSCANDLL_API GwyscanResult gwyscan_save_gwyddion_xyz(const double *xypos, const double *zdata, int n, int nchannels, int xres, int yres, const char **description, const char *filename);

GWYSCANDLL_API GwyscanResult gwyscan_add_gwyddion_xyz(const double *xypos, const double *zdata, int n, int nchannels, int xres, int yres, const char **description, const char *filename, const char *tmp_filename, int* n_total);
/**@}*/


/*!\ingroup GwyscanCreate
* @{
*/
GWYSCANDLL_API GwyscanResult gwyscan_generate_xyz_data(double *xypos, double *zdata, int n, double xreal, double yreal, double xoffset, double yoffset, int reserved);

GWYSCANDLL_API int gwyscan_create_path_regular(double *xypos, int xres, int yres, double xreal, double yreal, double xoffset, double yoffset, double angle, GwyscanDirection direction);

GWYSCANDLL_API int gwyscan_create_path_space_filling(double *xypos, int order, double xreal, double yreal, double xoffset, double yoffset);

GWYSCANDLL_API int gwyscan_create_path_space_filling_progressive(double *xypos, int order, int coarse_order, int skip_order, double xreal, double yreal, double xoffset, double yoffset);

GWYSCANDLL_API int gwyscan_create_path_spiral(double *xypos, int n, double xreal, double yreal, double xoffset, double yoffset);

GWYSCANDLL_API int gwyscan_create_path_random(double *xypos, int n, double xreal, double yreal, double xoffset, double yoffset);

GWYSCANDLL_API int gwyscan_create_path_octave_profiles(double *xypos, int order, int noctaves, int nprofiles, double xreal, double yreal, double xoffset, double yoffset);

GWYSCANDLL_API int gwyscan_create_path_octave_2d(double *xypos, int order, int noctaves, double xreal, double yreal, double xoffset, double yoffset);

GWYSCANDLL_API int gwyscan_create_path_wichmann_profiles(double *xypos, int n, int nprofiles, double xreal, double yreal, double xoffset, double yoffset);
/**@}*/


/*!\ingroup GwyscanModify
* @{
*/
GWYSCANDLL_API int gwyscan_modify_path_insert(double *xypos_dest, int *indexes, int *n_indexes, const double *xypos_source, int n_source, double xpos, double ypos, int interval);

GWYSCANDLL_API void gwyscan_modify_path_revert(double *xypos_dest, const double *xypos_source, int n);

GWYSCANDLL_API void gwyscan_modify_path_flip_horizontal(double *xypos_dest, const double *xypos_source, int n);

GWYSCANDLL_API void gwyscan_modify_path_flip_vertical(double *xypos_dest, const double *xypos_source, int n);
/**@}*/

// TODO: Refine Scan Path routines
#if(0)
/*!\ingroup GwyscanRefine
* @{
*/
GWYSCANDLL_API int gwyscan_refine_path_threshold(double* xypos, 
                                                 const double* xypos_prev, 
                                                 const double* zdata_prev, 
                                                 int n_prev, 
                                                 int channel_prev, 
                                                 const double* xypos_preprev,
                                                 const double* zdata_preprev, 
                                                 int n_preprev, 
                                                 int channel_preprev, 
                                                 double threshold, 
                                                 int xdiv, 
                                                 int ydiv, 
                                                 double xstep, 
                                                 double ystep, 
                                                 double xreal, 
                                                 double yreal, 
                                                 double xoffset, 
                                                 double yoffset);

GWYSCANDLL_API int gwyscan_refine_path_roughness(double* xypos, 
                                                 const double* xypos_prev, 
                                                 const double* zdata_prev, 
                                                 int n_prev, 
                                                 int channel_prev, 
                                                 double roughness, 
                                                 int xdiv, 
                                                 int ydiv, 
                                                 double xstep, 
                                                 double ystep, 
                                                 double xreal, 
                                                 double yreal, 
                                                 double xoffset, 
                                                 double yoffset);

GWYSCANDLL_API int gwyscan_refine_path_edges(double* xypos, 
                                             const double* xypos_prev, 
                                             const double* zdata_prev, 
                                             int n_prev, 
                                             int channel_prev, 
                                             double threshold, 
                                             double xstep, 
                                             double ystep, 
                                             double xreal, 
                                             double yreal, 
                                             double xoffset, 
                                             double yoffset);


/**@}*/
#endif
// TODO: Refine Scan Path routines


/*!\ingroup GwyscanInterpolation
* @{
*/
GWYSCANDLL_API GwyscanResult gwyscan_preview(const double *xypos, const double *zdata, int n, int channel, double *preview_image, double *preview_mask, int preview_width, int preview_height, double preview_xreal, double preview_yreal, double preview_xoffset, double preview_yoffset, double preview_angle);

GWYSCANDLL_API GwyscanResult gwyscan_preview_subset(const double *xypos, const double *zdata, int n, int from, int len, int channel_from, int channel_to, double **preview_images, double *preview_mask, int preview_width, int preview_height, double preview_xreal, double preview_yreal, double preview_xoffset, double preview_yoffset, double preview_angle);

GWYSCANDLL_API GwyscanResult gwyscan_preview_subset_all_channels(const double *xypos, const double *zdata, int n, int from, int len, int nchannels, double **preview_images, double *preview_mask, int preview_width, int preview_height, double preview_xreal, double preview_yreal, double preview_xoffset, double preview_yoffset, double preview_angle);

GWYSCANDLL_API GwyscanResult gwyscan_preview_average(const double *xypos, const double *zdata, int n, int channel, double *preview_image, double *preview_mask, int preview_width, int preview_height, double preview_xreal, double preview_yreal, double preview_xoffset, double preview_yoffset, double preview_angle);

GWYSCANDLL_API GwyscanResult gwyscan_preview_subset_average(const double *xypos, const double *zdata, int n, int from, int len, int channel_from, int channel_to, double **preview_images, double *preview_mask, int preview_width, int preview_height, double preview_xreal, double preview_yreal, double preview_xoffset, double preview_yoffset, double preview_angle);

GWYSCANDLL_API GwyscanResult gwyscan_preview_subset_average_all_channels(const double *xypos, const double *zdata, int n, int from, int len, int nchannels, double **preview_images, double *preview_mask, int preview_width, int preview_height, double preview_xreal, double preview_yreal, double preview_xoffset, double preview_yoffset, double preview_angle);
/**@}*/

/* vim: set cin et ts=4 sw=4 cino=>1s,e0,n0,f0,{0,}0,^0,\:1s,=0,g1s,h0,t0,+1s,c3,(0,u0 : */

