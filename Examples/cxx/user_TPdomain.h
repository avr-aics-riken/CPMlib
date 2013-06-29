#ifndef _TP_DOMAIN_BY_USER_
#define _TP_DOMAIN_BY_USER_

#include "cpm_TextParserDomain.h"

/** CPMの領域情報テキストパーサークラス */
class TPdomain : public cpm_TextParserDomain
{
public:
  /** コンストラクタ */
  TPdomain() : cpm_TextParserDomain()
  {
    // no operation
  }

  /** デストラクタ */
  virtual ~TPdomain()
  {
    // no operation
  }

  /** Read関数 */
  static cpm_GlobalDomainInfo* Read( std::string filename, int &errorcode );

  /** G_voxel,G_region,G_pitchからサイズを確定する */
  virtual int DecideGregion( bool bvox, int    vox[3]
                           , bool brgn, double rgn[3]
                           , bool bpch, double pch[3] );
};

#endif /* _TP_DOMAIN_BY_USER_ */

