#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cpm_ParaManager.h"
#include "user_TPdomain.h"

using namespace std;

//テストプログラムのメイン
int main( int argc, char **argv )
{
  int ret = 0;
  const char *ifname = "input.txt";
  cout << "input : " << ifname << endl;

  // 並列管理クラスのインスタンスと初期化
  cpm_ParaManager *paraMngr = cpm_ParaManager::get_instance(argc,argv);
  if( !paraMngr )
  {
    cerr << "PM instance : " << ret << endl;
    return CPM_ERROR_PM_INSTANCE;
  }

  // 領域分割情報の読み込み
#if 0
  cpm_GlobalDomainInfo *dInfo = cpm_TextParserDomain::Read(ifname,ret);
#else
  cpm_GlobalDomainInfo *dInfo = TPdomain::Read(ifname,ret);
#endif
  if( !dInfo )
  {
    cerr << "Read : pointer : " << ret << endl;
    return CPM_ERROR_INVALID_PTR;
  }
  if( ret != TP_NO_ERROR )
  {
    delete dInfo;
    cerr << "TextParser error : " << ret << endl;
    return ret;
  }

  if( paraMngr->GetMyRankID() == 0 )
  {
    cout << "vox=" << dInfo->GetVoxNum()[0] << ","
                   << dInfo->GetVoxNum()[1] << ","
                   << dInfo->GetVoxNum()[2] << endl;
    cout << "rgn=" << dInfo->GetRegion()[0] << ","
                   << dInfo->GetRegion()[1] << ","
                   << dInfo->GetRegion()[2] << endl;
    cout << "pch=" << dInfo->GetPitch()[0] << ","
                   << dInfo->GetPitch()[1] << ","
                   << dInfo->GetPitch()[2] << endl;
  }

  delete dInfo;
  return CPM_SUCCESS;
}

