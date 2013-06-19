/*
 * CPM - Cartesian Partition Manager
 *
 * Copyright (C) 2013 University of Tokyo.
 *
 */

/**
 * @file   cpm_EndianUtil.h
 * CPMエンディアンユーティリティヘッダーファイル
 * @author University of Tokyo
 * @date   2013/04/02
 */

#ifndef _CPM_ENDIAN_UTIL_H_
#define _CPM_ENDIAN_UTIL_H_

#include "cpm_Base.h"

/** CPMのエンディアンユーティリティ名前空間 */
namespace CPM_ENDIAN
{

  /** 16ビット(2バイト)変数のエンディアン変換
   * @param[inout] x 変換する変数
   */
  template<class X> CPM_INLINE void BSWAP16(X& x) {
    register unsigned char* _x_v = (unsigned char*)(&(x));
    unsigned char tmp;
    tmp = _x_v[0];
    _x_v[0] = _x_v[1];
    _x_v[1] = tmp;
  }

  /** 32ビット(4バイト)変数のエンディアン変換
   * @param[inout] x 変換する変数
   */
  template<class X> CPM_INLINE void BSWAP32(X& x) {
    register unsigned char* _x_v = (unsigned char*)(&(x));
    unsigned char tmp[2];
    tmp[0] = _x_v[0];
    tmp[1] = _x_v[1];
    _x_v[0] = _x_v[3];
    _x_v[1] = _x_v[2];
    _x_v[2] = tmp[1];
    _x_v[3] = tmp[0];
  }

  /** 64ビット(8バイト)変数のエンディアン変換
   * @param[inout] x 変換する変数
   */
  template<class X> CPM_INLINE void BSWAP64(X& x) {
    register unsigned char* _x_v = (unsigned char*)(&(x));
    unsigned char tmp[4];
    tmp[0] = _x_v[0];
    tmp[1] = _x_v[1];
    tmp[2] = _x_v[2];
    tmp[3] = _x_v[3];
    _x_v[0] = _x_v[7];
    _x_v[1] = _x_v[6];
    _x_v[2] = _x_v[5];
    _x_v[3] = _x_v[4];
    _x_v[4] = tmp[3];
    _x_v[5] = tmp[2];
    _x_v[6] = tmp[1];
    _x_v[7] = tmp[0];
  }

  /** 16ビット(2バイト)変数配列のエンディアン変換
   * @param[inout] a 変換する変数配列
   * @param[in]    n 配列要素数
   */
  template<class X, class Y> CPM_INLINE void SBSWAPVEC(X* a, Y n) {
    register unsigned int nn = (unsigned int)n;
    for(register unsigned int _i=0;_i<nn;_i++){
      register unsigned short _x_v = (unsigned short)a[_i];
      BSWAP16(_x_v);
      a[_i] = _x_v;
    }
  }


  /** 32ビット(4バイト)変数配列のエンディアン変換
   * @param[inout] a 変換する変数配列
   * @param[in]    n 配列要素数
   */
  template<class X, class Y> CPM_INLINE void BSWAPVEC(X* a, Y n) {
    register unsigned int nn = (unsigned int)n;
    for(register unsigned int _i=0;_i<nn;_i++){
      register unsigned int _x_v = (unsigned int)a[_i];
      BSWAP32(_x_v);
      a[_i] = _x_v;
    }
  }

  /** 64ビット(8バイト)変数配列のエンディアン変換
   * @param[inout] a 変換する変数配列
   * @param[in]    n 配列要素数
   */
  template<class X, class Y> CPM_INLINE void DBSWAPVEC(X* a, Y n) {
    register unsigned int nn = (unsigned int)n;
    for(register unsigned int _i=0;_i<nn;_i++){
      register unsigned long long _x_v = (unsigned long long)a[_i];
      BSWAP64(_x_v);
      a[_i] = _x_v;
    }
  }

  /** エンディアンチェックフラグ */
  enum EMatchType
  {
    UnKnown = 0 ///< 未定(フォーマット不明)
  , Match   = 1 ///< 一致
  , UnMatch = 2 ///< 不一致
  };

};

#endif /* _CPM_ENDIAN_UTIL_H_ */

