#include <chrono>
#include <cstring>
#include <iostream>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include <random>
#include <x86intrin.h>

#ifdef WIN32
#include <windows.h>

void ConvertTimeTToSystemTime(SYSTEMTIME* stime, time_t* time){
    tm my_tm;
    localtime_s(&my_tm, time);
    memset(stime, 0, sizeof(SYSTEMTIME));
    stime->wYear = my_tm.tm_year + 1900;
    stime->wMonth = my_tm.tm_mon + 1;
    stime->wDay = my_tm.tm_mday;
}

#define _cpuid __cpuid
#else
#include <cpuid.h>

int _cpuid(uint32_t* info, uint32_t func_id) {
    return __get_cpuid(func_id, info, info+1, info+2, info+3);
}
#endif

#define FILL_ZERO(NUM, ZEROS) std::setfill('0') << std::setw(ZEROS) << NUM
#define MAT_N 1024

#define MAT_MUL(I1, I2, I3) template <int N> \
    void mat_mul_##I1##I2##I3(int out[][N], int a[][N], int b[][N]) { \
        memset(out, 0, N * N * sizeof(int)); \
        for (int I1 = 0; I1 < N; I1++) { \
            for (int I2 = 0; I2 < N; I2++) \
                for (int I3 = 0; I3 < N; I3++) \
                    out[i][k] += a[i][k] * b[k][j]; \
        } \
    }

static int mat1024_a[MAT_N][MAT_N];
static int mat1024_b[MAT_N][MAT_N];
static int mat1024_out[MAT_N][MAT_N];


void fill_matrix1024() {
    std::uniform_int_distribution<int> randint(0, 4);
    std::default_random_engine eng;

    for (int i = 0; i < MAT_N; i++) {
        for (int j = 0; j < MAT_N; j++) {
            mat1024_a[i][j] = randint(eng);
            mat1024_b[i][j] = randint(eng);
        }
    }
}

void task_1() {
    time_t last_time = 0x7FFFFFFF;
    const auto l_time = gmtime(&last_time);
    std::cout << "Year: " << 1900 + l_time->tm_year
        << ", Month: " << FILL_ZERO(l_time->tm_mon + 1, 2)
        << ", Day: " << FILL_ZERO(l_time->tm_mday, 2) << std::endl;

    std::cout << "Full time: " << FILL_ZERO(l_time->tm_mday, 2)
        << "." << FILL_ZERO(l_time->tm_mon + 1, 2)
        << "." << 1900 + l_time->tm_year
        << " " << FILL_ZERO(l_time->tm_hour, 2)
        << ":" << FILL_ZERO(l_time->tm_min, 2)
        << ":" << FILL_ZERO(l_time->tm_sec, 2)
    << std::endl;

#ifdef WIN32
    SYSTEMTIME stime;
    ConvertTimeTToSystemTime(&stime, &last_time);
    printf("SYSTEMTIME: %4.4d.%2.2d.%2.2d\n", stime.wYear, stime.wMonth, stime.wDay);
#endif
}

time_t task_2_time() {
    const time_t start = time(nullptr);
    time_t finish = start;
    while (finish == start) {
        finish = time(nullptr);
    }

    return finish - start;
}

double task_2_omp() {
    const double start = omp_get_wtime();
    double finish = start;
    do {
        finish = omp_get_wtime();
    } while (finish - start == 0);

    return finish - start;
}

double task_2_chrono() {
    const std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
    std::chrono::high_resolution_clock::time_point finish = start;
    while ((finish - start).count() == 0) {
        finish = std::chrono::high_resolution_clock::now();
    }
    return std::chrono::duration_cast<std::chrono::duration<double>>(finish - start).count();
}

void task_2() {
    std::cout << "time: " << task_2_time() << " seconds" << std::endl;
    auto omp_res = task_2_omp();
    std::cout << "OpenMP: " << std::fixed << std::setprecision(10) << omp_res << " seconds (" << std::setprecision(0) << (omp_res * 1000000000) << " nanoseconds)" << std::endl;
    auto chrono_res = task_2_chrono();
    std::cout << "std::chrono: " << std::fixed << std::setprecision(10) << chrono_res << " seconds (" << std::setprecision(0) << (chrono_res * 1000000000) << " nanoseconds)" << std::endl;
}

MAT_MUL(i, j, k)
MAT_MUL(i, k, j)
MAT_MUL(j, k, i)
MAT_MUL(j, i, k)
MAT_MUL(k, i, j)
MAT_MUL(k, j, i)

double task_3_measure(int out[][MAT_N], int a[][MAT_N], int b[][MAT_N], void(*func)(int[][MAT_N], int[][MAT_N], int[][MAT_N])) {
    double time_min = MAXFLOAT;
    for (int i = 0; i < 5; i++){
        double time_cur = omp_get_wtime();
        func(out, a, b);
        time_cur = omp_get_wtime() - time_cur;
        if (time_cur < time_min)
            time_min = time_cur;
    }

    return time_min;
}

void task_3() {
    const double min_ijk = task_3_measure(mat1024_out, mat1024_a, mat1024_b, mat_mul_ijk<MAT_N>);
    printf("mul_ijk: c[0][0] = %d, c[N-1][N-1] = %d, time = %lg\n", mat1024_out[0][0], mat1024_out[MAT_N - 1][MAT_N - 1], min_ijk);

    const double min_ikj = task_3_measure(mat1024_out, mat1024_a, mat1024_b, mat_mul_ikj<MAT_N>);
    printf("mul_ikj: c[0][0] = %d, c[N-1][N-1] = %d, time = %lg\n", mat1024_out[0][0], mat1024_out[MAT_N - 1][MAT_N - 1], min_ikj);

    const double min_jki = task_3_measure(mat1024_out, mat1024_a, mat1024_b, mat_mul_jki<MAT_N>);
    printf("mul_jki: c[0][0] = %d, c[N-1][N-1] = %d, time = %lg\n", mat1024_out[0][0], mat1024_out[MAT_N - 1][MAT_N - 1], min_jki);

    const double min_jik = task_3_measure(mat1024_out, mat1024_a, mat1024_b, mat_mul_jik<MAT_N>);
    printf("mul_jik: c[0][0] = %d, c[N-1][N-1] = %d, time = %lg\n", mat1024_out[0][0], mat1024_out[MAT_N - 1][MAT_N - 1], min_jik);

    const double min_kij = task_3_measure(mat1024_out, mat1024_a, mat1024_b, mat_mul_kij<MAT_N>);
    printf("mul_kij: c[0][0] = %d, c[N-1][N-1] = %d, time = %lg\n", mat1024_out[0][0], mat1024_out[MAT_N - 1][MAT_N - 1], min_kij);

    const double min_kji = task_3_measure(mat1024_out, mat1024_a, mat1024_b, mat_mul_kji<MAT_N>);
    printf("mul_kji: c[0][0] = %d, c[N-1][N-1] = %d, time = %lg\n", mat1024_out[0][0], mat1024_out[MAT_N - 1][MAT_N - 1], min_kji);
}

void task_4() {
    uint64_t ll_min = 0x7FFFFFFFFFFFFFFF;
    uint32_t r[4];
    for (int i = 0; i < 5; i++) {
        _cpuid((uint32_t*)&r, 0);
        uint64_t ll_cur = __rdtsc();
        mat_mul_ikj(mat1024_out, mat1024_a, mat1024_b);
        _cpuid((uint32_t*)&r, 0);
        ll_cur = __rdtsc() - ll_cur;
        if (ll_cur < ll_min)
            ll_min = ll_cur;
    }

    printf("mul_ikj: c[0][0] = %d, c[N-1][N-1] = %d, time = %llu\n", mat1024_out[0][0], mat1024_out[MAT_N - 1][MAT_N - 1], ll_min);
}

void task_5() {
    double time_cur = 0;
    uint64_t ll_cur = 0;
    uint64_t ll_min = UINT64_MAX;
    uint32_t r[4];
    for (int i = 0; i < 500; i++) {
        _cpuid((uint32_t*)&r, 0);
        ll_cur = __rdtsc();
        time_cur += omp_get_wtime();
        _cpuid((uint32_t*)&r, 0);
        ll_cur = __rdtsc() - ll_cur;

        if (ll_cur < ll_min)
            ll_min = ll_cur;
    }
    printf("omp_get_wtime: ll_cur = %llu, time = %llu\n", ll_cur, ll_min);

    /*DWORD tick_min = 0xFFFFFFFF, tick_cur = 0;
    for (int i = 0; i < 500; i++) {
        __cpuid(r, 0);
        ll_cur = __rdtsc();
        tick_cur += GetTickCount();
        __cpuid(r, 0);
        ll_cur = __rdtsc() - ll_cur;
        if (ll_cur < ll_min)
            ll_min = ll_cur;
    }
    printf("GetTickCount : tick_cur = %d time = %I64d\n",
           tick_cur, ll_min);*/
}

void task_6() {
    constexpr int size = 100;

    int mat_a[size][size];
    int mat_b[size][size];
    int mat_out[size][size];

    int count = 0;
    const double start_time = omp_get_wtime();
    while (omp_get_wtime() - start_time < .2) {
        count++;
        mat_mul_ikj<size>(mat_out, mat_a, mat_b);
    }
    printf("mul_ijk: c[0][0] = %d, c[99][99] = %d, count = %d\n", mat_out[0][0], mat_out[99][99], count);
}

int main() {
    fill_matrix1024();

    task_1();
    std::cout << "\n";

    task_2();
    std::cout << "\n";

    task_3();
    std::cout << "\n";

    task_4();
    std::cout << "\n";

    task_5();
    std::cout << "\n";

    task_6();
    std::cout << "\n";

    return 0;
}
