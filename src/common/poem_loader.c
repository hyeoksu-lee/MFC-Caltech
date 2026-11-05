/* POEM host/device ABI loader. Device symbols are optional for CPU builds. */
#include <dlfcn.h>
#include <stdio.h>
#include <string.h>

static void *library;
static int (*initialize)(const char *);
static int (*constants)(double *);
static int (*state_dimension)(void);
static int (*decode)(int, const double *, double *, int *);
static int (*advance)(int, const double *, double, const double *, double *, double *, int *, double *, double *);
static int (*last_error)(char *, int);
static int (*device_initialize)(const char *, int);
static int (*device_decode)(int, const double *, double *, unsigned long long);
static int (*device_advance)(int, const int *, const double *, double, const double *,
                             double *, double *, double *, unsigned long long);
static char loader_error[2048];

static int load_library(const char *path) {
    if (!library) {
        library = dlopen(path, RTLD_NOW | RTLD_GLOBAL);
        if (!library) {
            snprintf(loader_error, sizeof(loader_error), "POEM dlopen: %s", dlerror());
            return 2;
        }
        initialize = (int (*)(const char *))dlsym(library, "poem_init");
        constants = (int (*)(double *))dlsym(library, "poem_constants");
        state_dimension = (int (*)(void))dlsym(library, "poem_state_dimension");
        decode = (int (*)(int,const double *,double *,int *))dlsym(library, "poem_decode");
        advance = (int (*)(int,const double *,double,const double *,double *,double *,int *,double *,double *))dlsym(library, "poem_step");
        last_error = (int (*)(char *,int))dlsym(library, "poem_last_error");
        device_initialize = (int (*)(const char *,int))dlsym(library, "poem_device_init");
        device_decode = (int (*)(int,const double *,double *,unsigned long long))dlsym(library, "poem_decode_device");
        device_advance = (int (*)(int,const int *,const double *,double,const double *,double *,double *,double *,unsigned long long))dlsym(library, "poem_step_device");
        if (!initialize || !state_dimension || !constants || !decode || !advance || !last_error) {
            snprintf(loader_error, sizeof(loader_error), "POEM backend ABI is incomplete");
            return 2;
        }
    }
    return 0;
}
int mfc_poem_init(const char *root, const char *path) {
    int code = load_library(path);
    return code ? code : initialize(root);
}
int mfc_poem_init_device(const char *root, const char *path, int device) {
    int code = load_library(path);
    if (code) return code;
    if (!device_initialize || !device_decode || !device_advance) {
        snprintf(loader_error, sizeof(loader_error), "POEM GPU build requires --backend jax-gpu");
        return 2;
    }
    return device_initialize(root, device);
}
int mfc_poem_decode_device(int n, const double *y, double *mom, unsigned long long stream) {
    return device_decode ? device_decode(n,y,mom,stream) : 2;
}
int mfc_poem_step_device(int n, const int *active, const double *p, double dt,
                        const double *y, double *next, double *ratio, double *root,
                        unsigned long long stream) {
    return device_advance ? device_advance(n,active,p,dt,y,next,ratio,root,stream) : 2;
}
int mfc_poem_state_dimension(void) { return state_dimension ? state_dimension() : 0; }
int mfc_poem_constants(double *v) { return constants(v); }
int mfc_poem_decode(int n, const double *y, double *mom, int *valid) {
    return decode(n, y, mom, valid);
}
int mfc_poem_step(int n, const double *p, double dt, const double *y,
                   double *next, double *mom, int *valid, double *ratio, double *root) {
    return advance(n, p, dt, y, next, mom, valid, ratio, root);
}
int mfc_poem_error(char *buffer, int n) {
    if (loader_error[0]) {
        snprintf(buffer, n, "%s", loader_error);
        return 0;
    }
    return last_error ? last_error(buffer, n) : 2;
}
