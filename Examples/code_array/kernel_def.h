#ifndef _KERNEL_DEF_H_
#define _KERNEL_DEF_H_
/*
 *  kernel_def.h
 *  kernel_test
 *
 *  Created by keno on 11/12/05.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#define VERS     5

#define JACOBI    1
#define SOR       2
#define SOR2SMA   3
#define PBICGSTAB 4
#define BICGSTAB  5
#define SOR2SMA_S 6
#define Vector_num 16

#define mark() printf("%s (%d):\n",__FILE__, __LINE__)


// FX10 profiler
#if defined __K_FPCOLL
#include "fjcoll.h"
#elif defined __FX_FAPP
#include "fj_tool/fjcoll.h"
#endif

// Performance Monitor
#include "PerfMonitor.h"

using namespace pm_lib;
using namespace std;

/**
 * @brief タイミング測定開始
 * @param [in] key ラベル
 */
inline void TIMING_start(PerfMonitor* pm, const string key)
{
    // PMlib Intrinsic profiler
    pm->start(key);
    
    const char* s_label = key.c_str();
    
    // Venus FX profiler
#if defined __K_FPCOLL
    start_collection( s_label );
#elif defined __FX_FAPP
    fapp_start( s_label, 0, 0);
#endif
}

/**
 * @brief タイミング測定終了
 * @param [in] key             ラベル
 * @param [in] flopPerTask    「タスク」あたりの計算量/通信量(バイト) (ディフォルト0)
 * @param [in] iterationCount  実行「タスク」数 (ディフォルト1)
 */
inline void TIMING_stop(PerfMonitor* pm, const string key, double flopPerTask=0.0, int iterationCount=1)
{
    // Venus FX profiler
    const char* s_label = key.c_str();
    
#if defined __K_FPCOLL
    stop_collection( s_label );
#elif defined __FX_FAPP
    fapp_stop( s_label, 0, 0);
#endif
    
    // PMlib Intrinsic profiler
    pm->stop(key, flopPerTask, (unsigned)iterationCount);
}
#endif
