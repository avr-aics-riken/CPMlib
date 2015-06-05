#include "cpm_Base.h"
#include "cpm_ParaManager.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mconvp_LMR_main.h"

using namespace std;

extern "C"
{
  void
  diffusion_( int *ischeme, int *stpMax, double *dt, double *cf, int *itrMax, double *eps, double *omg, int *ret );
}

#define _FLUSH_PG(_S,_PG) \
{ \
  paraMngr->Barrier(_PG); \
  _S << flush; \
  paraMngr->Barrier(_PG); \
}

#define _FLUSH(_S) _FLUSH_PG(_S,0)


int main( int argc, char **argv )
{
  int ret = 0;

  // 並列管理クラスのインスタンスと初期化
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance(argc,argv,CPM_DOMAIN_LMR);
  if( !paraMngr ) return CPM_ERROR_PM_INSTANCE;

  // 引数のチェック
  if( argc != 2 )
  {
    cout << "argument error" << endl;
    return 9;
  }

  // 時間計測開始
  double ts = cpm_Base::GetWTime();

  if( paraMngr->GetMyRankID() == 0 )
  {
    cout << "CPMlib Version " <<
    cpm_Base::getVersionInfo() << endl;
    cout << "CPMlib Revision " <<
    cpm_Base::getRevisionInfo() << endl;
  }

  // 入力ファイル
  const char *ifname = argv[1];
  cout << "input file : " << ifname << endl;
  InputParam input(ifname, ret);
  if( ret != 0 )
  {
    cout << "input file read error" << endl;
    return ret;
  }

  // 領域分割
  if( (ret = paraMngr->VoxelInit_LMR( ifname, 2, 1 )) != CPM_SUCCESS )
  {
    cerr << "VoxelInit_LMR error [" << ifname << "] : " << ret << endl;
    return ret;
  }

  // 自ランクの情報を取得
  int myrank = paraMngr->GetMyRankID();
  int nrank  = paraMngr->GetNumRank();
  const int    *vox = paraMngr->GetLocalVoxelSize();
  const double *pch = paraMngr->GetPitch();
/*
  _FLUSH(cout);
  for( int i=0;i<nrank;i++ )
  {
    if( i==myrank )
    {
      cout << "[" << myrank << "/" << nrank << "] voxel size  = " << vox[0] << "," << vox[1] << "," << vox[2] << endl
           << "[" << myrank << "/" << nrank << "] voxel pitch = " << pch[0] << "," << pch[1] << "," << pch[2] << endl;
    }
    _FLUSH(cout);
  }
*/

  // parameter
  string scheme = input.GetString("/MCONV/Scheme");
  int    stpMax = input.GetInt   ("/MCONV/stpMax");
  double dt     = input.GetDouble("/MCONV/dt");
  double cf     = input.GetDouble("/MCONV/cf");
  int itrMax = 0;
  double eps = 0.0;
  double omg = 0.0;
  int ischeme = 1;
  if( InputParam::strCompare(scheme, "jacobi") == 0 )
  {
    ischeme = 2;
    itrMax = input.GetInt   ("/MCONV/itrMax");
    eps    = input.GetDouble("/MCONV/eps");
    omg    = input.GetDouble("/MCONV/omg");
  }
/*
  if( myrank == 0 )
  {
    cout << "tp file : " << tpFile << endl;
    cout << "scheme  : " << scheme << endl;
    cout << "stpMax  : " << stpMax << endl;
    cout << "dt      : " << dt << endl;
    cout << "cf      : " << cf << endl;
    if( ischeme == 2 )
    {
      cout << "itrMax  : " << itrMax << endl;
      cout << "eps     : " << eps << endl;
      cout << "omg     : " << omg << endl;
    } 
  }
*/
  _FLUSH(cout);

  // diffusion
  diffusion_( &ischeme, &stpMax, &dt, &cf, &itrMax, &eps, &omg, &ret );

  return ret;
}


