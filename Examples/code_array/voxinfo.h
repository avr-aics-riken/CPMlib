/*
 *  voxinfo.h
 *  kernel_test
 *
 *  Created by keno on 11/12/05.
 *  Copyright 2011 __iis__. All rights reserved.
 *
 */

#ifndef _VOXINFO_H_
#define _VOXINFO_H_

#include <string>
#include "FBdefine.h"

class VoxInfo {
  public:
  VoxInfo() {}
  ~VoxInfo() {}
  
protected:
  unsigned encAS         (int* m_sz, int gc, int* pad_size, unsigned* bx);
  
  void encPbit           (int* m_sz, int gc, int* pad_size, unsigned* bx);
  void encPbit_OBC       (int* m_sz, int gc, int* pad_size, int face, unsigned* bx, std::string key, bool dir);
  
  /**
   @fn inline unsigned offBit(unsigned idx, const unsigned shift)
   @brief BCindexの第shiftビットをOFFにする
   @retval エンコードした値
   @param idx BCindex
   @param shift 指定ビット
   */
  inline unsigned offBit(unsigned idx, const unsigned shift) {
    return ( idx & (~(0x1<<shift)) );
  }
  
  /**
   @fn inline unsigned onBit(unsigned idx, const unsigned shift)
   @brief BCindexの第shiftビットをONにする
   @retval エンコードした値
   @param idx BCindex
   @param shift 指定ビット
   */
  inline unsigned onBit(unsigned idx, const unsigned shift) {
    return ( idx | (0x1<<shift) );
  }
  
  public:

  void setBCIndexP       (int* m_sz, int gc, int* pad_size, unsigned* bcp);
  
};


#endif // _VOXINFO_H_
