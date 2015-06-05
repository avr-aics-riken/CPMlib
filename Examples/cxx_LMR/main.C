/*
 * CPMlib - Cartesian Partition Manager Library
 *
 * Copyright (C) 2012-2014 Institute of Industrial Science, The University of Tokyo.
 * All rights reserved.
 *
 * Copyright (c) 2014-2015 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#include "cpm_ParaManager.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include "cpm_TextParserDomain.h"

using namespace std;

#define _FLUSH_PG(_S,_PG) \
{ \
  paraMngr->Barrier(_PG); \
  _S << flush; \
  paraMngr->Barrier(_PG); \
}

#define _FLUSH(_S) _FLUSH_PG(_S,0)

//テストプログラムのメイン
int main( int argc, char **argv )
{
  int ret = 0;

  // 並列管理クラスのインスタンスと初期化
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance(argc,argv,CPM_DOMAIN_LMR);
  if( !paraMngr ) return CPM_ERROR_PM_INSTANCE;

  // 時間計測開始
  double ts = cpm_Base::GetWTime();

  if( paraMngr->GetMyRankID() == 0 )
  {
    cout << "CPMlib Version " << 
    cpm_Base::getVersionInfo() << endl;
  }

  // 入力ファイルリスト
  vector<const char*> ifname;
  for( int i=1;i<argc;i++ )
  {
    ifname.push_back(argv[i]);
  }

  // ちょっとした情報のプリント
  _FLUSH(cout);
  if( paraMngr->GetMyRankID()==0 )
  {
    //入力となる領域分割ファイル
    cout << "number of input files=" << ifname.size() << endl;
    for( int i=0;i<ifname.size();i++ )
    {
      cout << "  " << i << " : " << ifname[i] << endl;
    }
  }
  _FLUSH(cout);

  // 領域分割
  for( int i=0;i<ifname.size();i++ )
  {
    int procGrpNo = 0;
    if( i>0 )
    {
      int nproc = cpm_ParaManagerLMR::GetNumLeaf( ifname[i] );
      int *proclist = new int[nproc];
      for( int j=0;j<nproc;j++ ) proclist[j] = j;
      procGrpNo = paraMngr->CreateProcessGroup( nproc, proclist );
    }

    // 領域分割
    if( (ret = paraMngr->VoxelInit_LMR( ifname[i], 3, 6, procGrpNo )) != CPM_SUCCESS )
    {
      cerr << "VoxelInit_LMR error [" << ifname[i] << "] : " << ret << endl;
      return ret;
    }
  }

#ifdef _DEBUG
  paraMngr->printVoxelInfo();
  paraMngr->printVoxelInfo(paraMngr->GetMyRankID());
  _FLUSH(cout);
#endif

  ////// MPI test //////
  #include "include/_mpi_test.h"

  paraMngr->Barrier();
  double elapse = cpm_Base::GetWSpanTime(ts);
  if( paraMngr->GetMyRankID()==0 ) cout << endl << "### elapse = "<< elapse << " [sec]" << endl;

  return CPM_SUCCESS;
}

