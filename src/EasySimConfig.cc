
/*
  EasySimConfig.cc

  implementation file for class EasySimConfig

*/


#include "EasySimConfig.h"

//----------------------------------------------------------------
/*
  singleton bouzin
*/
//----------------------------------------------------------------
EasySimConfig::EasySimConfig()
{

}
EasySimConfig *EasySimConfig::_instance = new EasySimConfig();
EasySimConfig *EasySimConfig::Instance()
{
  return _instance;
}
EasySimConfig* theConfig()
{
  return EasySimConfig::Instance();
}
