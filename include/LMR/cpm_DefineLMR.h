/*
###################################################################################
#
# CPMlib - Computational space Partitioning Management library
#
# Copyright (c) 2012-2014 Institute of Industrial Science (IIS), The University of Tokyo.
# All rights reserved.
#
# Copyright (c) 2014-2016 Advanced Institute for Computational Science (AICS), RIKEN.
# All rights reserved.
#
# Copyright (c) 2016-2017 Research Institute for Information Technology (RIIT), Kyushu University.
# All rights reserved.
#
###################################################################################
 */

/**
 * @file   cpm_DefineLMR.h
 * CPM-LMR用の定義マクロ記述ヘッダーファイル
 * @date   2012/05/31
 */
#ifndef _CPM_DEFINE_LMR_H_
#define _CPM_DEFINE_LMR_H_

#include "cpm_Define.h"


/** 3次元インデクス(i,j,k,iLeaf) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _IL ローカルリーフ順番号(0～)
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_S3D_LMR(_I,_J,_K,_IL,_NI,_NJ,_NK,_VC) (_IDX_S4D(_I,_J,_K,_IL,_NI,_NJ,_NK,_VC))

/** 4次元インデクス(i,j,k,n,iLeaf) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _N  成分インデクス
 *  @param[in] _IL ローカルリーフ順番号(0～)
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _NN 成分数
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_S4D_LMR(_I,_J,_K,_N,_IL,_NI,_NJ,_NK,_NN,_VC) \
( (long long)(_NN) * (long long)(_NI+2*(_VC)) * (long long)(_NJ+2*(_VC)) * (long long)(_NK+2*(_VC)) \
* (long long)(_IL) \
+ _IDX_S4D(_I,_J,_K,_N,_NI,_NJ,_NK,_VC) \
)

/** 3次元インデクス(i,j,k,3,iLeaf) -> 1次元インデクス変換マクロ
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _N  成分インデクス
 *  @param[in] _IL ローカルリーフ順番号(0～)
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 */
#define _IDX_V3D_LMR(_I,_J,_K,_N,_IL,_NI,_NJ,_NK,_VC) (_IDX_S4D_LMR(_I,_J,_K,_N,_IL,_NI,_NJ,_NK,3,_VC))

/** 4次元インデクス(n,i,j,k,iLeaf) -> 1次元インデクス変換マクロ
 *  @param[in] _N  成分インデクス
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _IL ローカルリーフ順番号(0～)
 *  @param[in] _NN 成分数
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_S4DEX_LMR(_N,_I,_J,_K,_IL,_NN,_NI,_NJ,_NK,_VC) \
( (long long)(_NN) * (long long)(_NI+2*(_VC)) * (long long)(_NJ+2*(_VC)) * (long long)(_NK+2*(_VC)) \
* (long long)(_IL) \
+ _IDX_S4DEX(_N,_I,_J,_K,_NN,_NI,_NJ,_NK,_VC) \
)

/** 3次元インデクス(3,i,j,k,iLeaf) -> 1次元インデクス変換マクロ
 *  @param[in] _N  成分インデクス
 *  @param[in] _I  i方向インデクス
 *  @param[in] _J  j方向インデクス
 *  @param[in] _K  k方向インデクス
 *  @param[in] _IL ローカルリーフ順番号(0～)
 *  @param[in] _NI i方向インデクスサイズ
 *  @param[in] _NJ j方向インデクスサイズ
 *  @param[in] _NK k方向インデクスサイズ
 *  @param[in] _VC 仮想セル数
 *  @return 1次元インデクス
 */
#define _IDX_V3DEX_LMR(_N,_I,_J,_K,_IL,_NI,_NJ,_NK,_VC) (_IDX_S4DEX_LMR(_N,_I,_J,_K,_IL,3,_NI,_NJ,_NK,_VC))





#endif /* _CPM_DEFINE_LMR_H_ */
