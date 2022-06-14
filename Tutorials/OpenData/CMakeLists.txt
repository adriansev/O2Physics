# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.


o2physics_add_dpl_workflow(flow-analysis
                  SOURCES flowAnalysis.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)


target_include_directories(treecalibv0zdc PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
# Then generate dictionaries and add them as a dependency of the executable (via the MODULE parameter):
ROOT_GENERATE_DICTIONARY(TreeCalibV0ZDC_dict TreeCalibV0ZDC.h MODULE TreeCalibV0ZDC LINKDEF TreeCalibV0ZDC_linkdef.h)

o2physics_add_dpl_workflow(treecalibv0zdc
                  SOURCES TreeCalibV0ZDC.cxx
                  PUBLIC_LINK_LIBRARIES O2::Framework O2Physics::AnalysisCore
                  COMPONENT_NAME AnalysisTutorial)

