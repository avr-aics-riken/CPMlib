#ifndef _CPM_TESTPROG_DEBUG_H_
#define _CPM_TESTPROG_DEBUG_H_

#include "cpm_ParaManager.h"
#include "cpm_TextParserDomain.h"

//読み込んだ領域情報のデバッグライト
void printDomainInfo( cpm_GlobalDomainInfo* dInfo, cpm_ParaManager *pMngr, int procGrp=0 )
{
  const double *org = dInfo->GetOrigin();
  const double *pch = dInfo->GetPitch();
  const double *rgn = dInfo->GetRegion();
  const int    *vox = dInfo->GetVoxNum();
  const int    *div = dInfo->GetDivNum();

  const int    rank = pMngr->GetMyRankID(procGrp);
  const int    *pos = pMngr->GetDivPos(procGrp);

  std::cout << "####### read parameters ########" << std::endl;
  std::cout << " G_org      = " << org[0] << "," << org[1] << "," << org[2] << std::endl;
  std::cout << " G_voxel    = " << vox[0] << "," << vox[1] << "," << vox[2] << std::endl;
  std::cout << " G_pitch    = " << pch[0] << "," << pch[1] << "," << pch[2] << std::endl;
  std::cout << " G_region   = " << rgn[0] << "," << rgn[1] << "," << rgn[2] << std::endl;
  std::cout << " G_div      = " << div[0] << "," << div[1] << "," << div[2] << std::endl;
  std::cout << " #subdomain = " << dInfo->GetSubdomainNum() << std::endl;
  for( size_t i=0;i<dInfo->GetSubdomainArraySize();i++ )
  {
    const cpm_ActiveSubdomainInfo *dom = dInfo->GetSubdomainInfo(i);
    const int *pos  = dom->GetPos();
    std::cout << "  domain" << i
              << "  pos="   << pos[0]  << "," << pos[1]  << "," << pos[2]
              << std::endl;
  }
  std::cout << " procGrp = " << procGrp << std::endl;
  std::cout << " my rank = " << rank << std::endl;
  std::cout << " my pos  = " << pos[0] << "," << pos[1] << "," << pos[2] << std::endl;
  std::cout << std::endl;

}

#endif /* _CPM_TESTPROG_DEBUG_H_ */
