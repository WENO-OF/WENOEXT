/*---------------------------------------------------------------------------* \
License
    This file is part of libWallModelledLES.

    libWallModelledLES is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    libWallModelledLES is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with libWallModelledLES.
    If not, see <http://www.gnu.org/licenses/>.

Class


Description
    Based on the OpenFOAM version, generates macros that can be used to make
    the code support multiple versions via #ifdef statements.
    The file is adopted with minor changes from swak4foam, see list of
    contributors below.

SourceFiles

Contributors/Copyright:
    2012-2018 Bernhard F.W. Gschaider <bgschaid@hfd-research.com>
    2015-2016 Bruno Santos <wyldckat@gmail.com>
    2016-2018 Mark Olesen <Mark.Olesen@esi-group.com>
    2019      Timofey Mukha

\*---------------------------------------------------------------------------*/

#ifndef WenoMacroHeader_H
#define WenoMacroHeader_H

#include "foamVersion4weno.H"

#define VERSION_NR(major,minor,patch) (     \
    10000 * major + \
    100   * minor + \
    patch )

#define VERSION_NR2(major,minor) (     \
    10000 * major + \
    100   * minor + \
    99 )

#define FOAM_VERSION4WENO VERSION_NR(FOAM_VERSION4WENO_MAJOR,FOAM_VERSION4WENO_MINOR,FOAM_VERSION4WENO_PATCH_NUM)

#if __cplusplus > 199711L
#define FOAM_CXX_COMPILER_ALLOWS_NON_POD_IN_UNION
#endif

#if __cplusplus < 201103L
#define FOAM_CXX_COMPILER_HAS_NO_NULLPTR
#ifndef nullptr
#define nullptr NULL
#endif
#endif

#if (OPENFOAM_COM >= 1606) \
    || (defined(FOAM_VERSION4WENO_IS_ORG) && FOAM_VERSION4WENO>VERSION_NR2(3,0))
#define FOAM_NEW_TMP_RULES
#endif

#if (OPENFOAM_COM >= 1612) \
    || (defined(FOAM_VERSION4WENO_IS_ORG) && FOAM_VERSION4WENO>VERSION_NR2(3,0))
#define FOAM_NEW_GEOMFIELD_RULES
#endif

#if (OPENFOAM_COM >= 1912) \
    || (defined(FOAM_VERSION4WENO_IS_ORG) && FOAM_VERSION4WENO>VERSION_NR2(6,0))
#define FOAM_MOVE_CONSTRUCTOR
#endif

#if (OPENFOAM_COM >= 1706) \
    || (defined(FOAM_VERSION4WENO_IS_ORG) && FOAM_VERSION4WENO>=VERSION_NR(5,0,0))
#define FOAM_PSTREAM_COMMSTYPE_IS_ENUMCLASS
#endif

#if (OPENFOAM_COM >= 1706) \
    || (defined(FOAM_VERSION4WENO_IS_ORG) && FOAM_VERSION4WENO>=VERSION_NR(5,0,0))
#define FOAM_PSTREAM_COMMSTYPE_IS_ENUMCLASS
#endif

#if ((defined(FOAM_VERSION4WENO_IS_ORG) && FOAM_VERSION4WENO>=VERSION_NR(8,0,0)))
#define FOAM_NEW_COUPLED_PATCHES
#endif



#endif

// ************************************************************************* //
