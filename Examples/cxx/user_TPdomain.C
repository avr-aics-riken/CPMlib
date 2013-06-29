#include "user_TPdomain.h"

using namespace std;

// Read関数
cpm_GlobalDomainInfo*
TPdomain::Read( std::string filename, int &errorcode )
{
  // インスタンス
  TPdomain tp;

  // 読み込みメイン
  return tp.ReadMain( filename, errorcode );
}

// G_voxel,G_region,G_pitchからサイズを確定する
int
TPdomain::DecideGregion( bool bvox, int    vox[3]
                       , bool brgn, double rgn[3]
                       , bool bpch, double pch[3] )
{
  cout << __PRETTY_FUNCTION__ << endl;
  // G_voxelが記述されている
  if( bvox )
  {
    // G_region優先
    if( brgn )
    {
      if( rgn[0] <= 0.0 || rgn[0] <= 0.0 || rgn[0] <= 0.0 )
      {
        return CPM_ERROR_TP_INVALID_G_RGN;
      }
      for( int i=0;i<3;i++ ) pch[i] = rgn[i] / double(vox[i]);
    }
    // G_regionの記述が無く、G_pitchが記述されている
    else if( bpch )
    {
      if( pch[0] <= 0.0 || pch[0] <= 0.0 || pch[0] <= 0.0 )
      {
        return CPM_ERROR_TP_INVALID_G_PITCH;
      }
      for( int i=0;i<3;i++ ) rgn[i] = pch[i] * double(vox[i]);
    }
    // G_regionもG_pitchも記述無し
    else
    {
      return CPM_ERROR_TP_INVALID_G_RGN;
    }
  }
  // G_voxelが記述されていない
  else
  {
    if( !brgn || !bpch )
    {
      return CPM_ERROR_TP_INVALID_G_RGN;
    }
    // 切り上げ
    for( int i=0;i<3;i++ )
    {
      vox[i] = (int)ceil(rgn[i]/pch[i]);
    }
  }

  return CPM_SUCCESS;
}

