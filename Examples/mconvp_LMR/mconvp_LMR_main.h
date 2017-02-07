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

#pragma once

#include "cpm_TextParser.h"
#include "cpm_PathUtil.h"

class InputParam : cpm_Base
{
private:
  const char *m_ifname;
  TextParser *m_tp;
  InputParam()
  {
    m_ifname = NULL;
    m_tp = NULL;
  }

public:
  InputParam(const char *ifname, int &ret)
  {
    ret = 0;

    m_ifname = ifname;

    m_tp = new TextParser();
    if( !m_tp )
    {
      ret = TP_ERROR;
      return;
    }

    if( (ret = m_tp->remove()) != TP_NO_ERROR ) return;

    if( (ret = m_tp->read(m_ifname)) != TP_NO_ERROR ) return;
  }

  ~InputParam()
  {
    if( m_tp ) delete m_tp;
  }

  std::string GetString( std::string label )
  {
    std::string value;
    if( m_tp->getValue(label, value) != TP_NO_ERROR ) return std::string("");
    return value;
  }

  int GetInt( std::string label )
  {
    int ret = 0;
    std::string value;
    if( m_tp->getValue(label, value) != TP_NO_ERROR ) return 0;
    return m_tp->convertInt(value, &ret);
  }

  double GetDouble( std::string label )
  {
    int ret = 0;
    std::string value;
    if( m_tp->getValue(label, value) != TP_NO_ERROR ) return 0.0;
    return m_tp->convertDouble(value, &ret);
  }

  static int strCompare( std::string str1, std::string str2, bool ignorecase=true )
  {
    InputParam base;
    return base.cpm_strCompare( str1, str2, ignorecase );
  }

};
