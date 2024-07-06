/* -----------------------------------------------------------------------------
Software License for The Fraunhofer FDK AAC Codec Library for Android

© Copyright  1995 - 2018 Fraunhofer-Gesellschaft zur Förderung der angewandten
Forschung e.V. All rights reserved.

 1.    INTRODUCTION
The Fraunhofer FDK AAC Codec Library for Android ("FDK AAC Codec") is software
that implements the MPEG Advanced Audio Coding ("AAC") encoding and decoding
scheme for digital audio. This FDK AAC Codec software is intended to be used on
a wide variety of Android devices.

AAC's HE-AAC and HE-AAC v2 versions are regarded as today's most efficient
general perceptual audio codecs. AAC-ELD is considered the best-performing
full-bandwidth communications codec by independent studies and is widely
deployed. AAC has been standardized by ISO and IEC as part of the MPEG
specifications.

Patent licenses for necessary patent claims for the FDK AAC Codec (including
those of Fraunhofer) may be obtained through Via Licensing
(www.vialicensing.com) or through the respective patent owners individually for
the purpose of encoding or decoding bit streams in products that are compliant
with the ISO/IEC MPEG audio standards. Please note that most manufacturers of
Android devices already license these patent claims through Via Licensing or
directly from the patent owners, and therefore FDK AAC Codec software may
already be covered under those patent licenses when it is used for those
licensed purposes only.

Commercially-licensed AAC software libraries, including floating-point versions
with enhanced sound quality, are also available from Fraunhofer. Users are
encouraged to check the Fraunhofer website for additional applications
information and documentation.

2.    COPYRIGHT LICENSE

Redistribution and use in source and binary forms, with or without modification,
are permitted without payment of copyright license fees provided that you
satisfy the following conditions:

You must retain the complete text of this software license in redistributions of
the FDK AAC Codec or your modifications thereto in source code form.

You must retain the complete text of this software license in the documentation
and/or other materials provided with redistributions of the FDK AAC Codec or
your modifications thereto in binary form. You must make available free of
charge copies of the complete source code of the FDK AAC Codec and your
modifications thereto to recipients of copies in binary form.

The name of Fraunhofer may not be used to endorse or promote products derived
from this library without prior written permission.

You may not charge copyright license fees for anyone to use, copy or distribute
the FDK AAC Codec software or your modifications thereto.

Your modified versions of the FDK AAC Codec must carry prominent notices stating
that you changed the software and the date of any change. For modified versions
of the FDK AAC Codec, the term "Fraunhofer FDK AAC Codec Library for Android"
must be replaced by the term "Third-Party Modified Version of the Fraunhofer FDK
AAC Codec Library for Android."

3.    NO PATENT LICENSE

NO EXPRESS OR IMPLIED LICENSES TO ANY PATENT CLAIMS, including without
limitation the patents of Fraunhofer, ARE GRANTED BY THIS SOFTWARE LICENSE.
Fraunhofer provides no warranty of patent non-infringement with respect to this
software.

You may use this FDK AAC Codec software or modifications thereto only for
purposes that are authorized by appropriate patent licenses.

4.    DISCLAIMER

This FDK AAC Codec software is provided by Fraunhofer on behalf of the copyright
holders and contributors "AS IS" and WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES,
including but not limited to the implied warranties of merchantability and
fitness for a particular purpose. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE for any direct, indirect, incidental, special, exemplary,
or consequential damages, including but not limited to procurement of substitute
goods or services; loss of use, data, or profits, or business interruption,
however caused and on any theory of liability, whether in contract, strict
liability, or tort (including negligence), arising in any way out of the use of
this software, even if advised of the possibility of such damage.

5.    CONTACT INFORMATION

Fraunhofer Institute for Integrated Circuits IIS
Attention: Audio and Multimedia Departments - FDK AAC LL
Am Wolfsmantel 33
91058 Erlangen, Germany

www.iis.fraunhofer.de/amm
amm-info@iis.fraunhofer.de
----------------------------------------------------------------------------- */

/******************* Library for basic calculation routines ********************

   Author(s):   Dangwi Xu @ sophgo

   Description: mdct RISC-V assembler replacements.

*******************************************************************************/

#if defined(__riscv_xthead) && __riscv_v == 7000
#include <riscv_vector.h>

#ifndef FUNCTION_mdct_block
#define FUNCTION_mdct_block

INT mdct_block(H_MDCT hMdct, const INT_PCM *RESTRICT timeData,
               const INT noInSamples, FIXP_DBL *RESTRICT mdctData,
               const INT nSpec, const INT tl, const FIXP_WTP *pRightWindowPart,
               const INT fr, SHORT *pMdctData_e) {
  int i, n;
  /* tl: transform length
     fl: left window slope length
     nl: left window slope offset
     fr: right window slope length
     nr: right window slope offset
     See FDK_tools/doc/intern/mdct.tex for more detail. */
  int fl, nl, nr;
  const FIXP_WTP *wls, *wrs;

  wrs = pRightWindowPart;

  /* Detect FRprevious / FL mismatches and override parameters accordingly */
  if (hMdct->prev_fr ==
      0) { /* At start just initialize and pass parameters as they are */
    hMdct->prev_fr = fr;
    hMdct->prev_wrs = wrs;
    hMdct->prev_tl = tl;
  }

  /* Derive NR */
  nr = (tl - fr) >> 1;

  /* Skip input samples if tl is smaller than block size */
  timeData += (noInSamples - tl) >> 1;

  /* windowing */
  for (n = 0; n < nSpec; n++) {
    INT mdctData_e = 1 + 1;

    /* Derive left parameters */
    wls = hMdct->prev_wrs;
    fl = hMdct->prev_fr;
    nl = (tl - fl) >> 1;

    for (i = 0; i < nl; i++) {
#if SAMPLE_BITS == DFRACT_BITS /* SPC_BITS and DFRACT_BITS should be equal. */
      mdctData[(tl / 2) + i] = -((FIXP_DBL)timeData[tl - i - 1] >> (1));
#else
      mdctData[(tl / 2) + i] = -(FIXP_DBL)timeData[tl - i - 1]
                               << (DFRACT_BITS - SAMPLE_BITS - 1); /* 0(A)-Br */
#endif
    }

    {
      INT vl, vlmax;
      vuint16m1_t vtimeidx, vtimelen;
      vint16m1_t  vtimelo, vtimehi, vim, vre;
      vint32m2_t  vmdct, vtmp0, vtmp1;

      vlmax = vsetvlmax_e16m1();

      vtimelen = vmv_v_x_u16m1((tl - nl - 1) * 2, vlmax);
      vtimeidx = vid_v_u16m1(vlmax);
      vtimeidx = vsll_vx_u16m1(vtimeidx, 1, vlmax);
      vtimeidx = vsub_vv_u16m1(vtimelen, vtimeidx, vlmax);

      for (i = 0; i < fl / 2; i += vl) {
        vl = vsetvl_e16m1(fl / 2 - i);

        vtimelo = vle16_v_i16m1(&timeData[i + nl], vl);
        vtimehi = vlxh_v_i16m1(timeData, vtimeidx, vl);
        vre = vlse16_v_i16m1(((int16_t*)&wls[i]) + 0, sizeof(FIXP_WTP), vl);
        vim = vlse16_v_i16m1(((int16_t*)&wls[i]) + 1, sizeof(FIXP_WTP), vl);

        vtmp0 = vwmul_vv_i32m2(vtimelo, vim, vl);
        vtmp1 = vwmul_vv_i32m2(vtimehi, vre, vl);

        vmdct = vsub_vv_i32m2(vtmp0, vtmp1, vl);
        vtimeidx = vsub_vx_u16m1(vtimeidx, vl * 2, vl);

        vse32_v_i32m2(&mdctData[(tl / 2) + i + nl], vmdct, vl);
      }
    }

    for (i = 0; i < nr; i++) {
#if SAMPLE_BITS == \
    DFRACT_BITS /* This should be SPC_BITS instead of DFRACT_BITS. */
      mdctData[(tl / 2) - 1 - i] = -((FIXP_DBL)timeData[tl + i] >> (1));
#else
      mdctData[(tl / 2) - 1 - i] =
          -(FIXP_DBL)timeData[tl + i]
          << (DFRACT_BITS - SAMPLE_BITS - 1); /* -C flipped at placing */
#endif
    }

    {
      INT vl, vlmax;
      vuint16m1_t vtimeidx, vtimelen;
      vuint32m2_t vmdctidx, vmdctlen;
      vint16m1_t  vtimelo, vtimehi, vim, vre;
      vint32m2_t  vmdct, vtmp0, vtmp1;

      vlmax = vsetvlmax_e16m1();

      vtimelen = vmv_v_x_u16m1(((tl * 2) - nr - 1) * 2, vlmax);
      vtimeidx = vid_v_u16m1(vlmax);
      vtimeidx = vsll_vx_u16m1(vtimeidx, 1, vlmax);
      vtimeidx = vsub_vv_u16m1(vtimelen, vtimeidx, vlmax);

      vmdctlen = vmv_v_x_u32m2(((tl / 2) - nr - 1) * 4, vlmax);
      vmdctidx = vid_v_u32m2(vlmax);
      vmdctidx = vsll_vx_u32m2(vmdctidx, 2, vlmax);
      vmdctidx = vsub_vv_u32m2(vmdctlen, vmdctidx, vlmax);

      for (i = 0; i < fr / 2; i += vl) {
        vl = vsetvl_e16m1(fr / 2 - i);

        vtimelo = vle16_v_i16m1(&timeData[tl + nr + i], vl);
        vtimehi = vlxh_v_i16m1(timeData, vtimeidx, vl);
        vre  = vlse16_v_i16m1(((int16_t*)&wrs[i]) + 0, sizeof(FIXP_WTP), vl);
        vim  = vlse16_v_i16m1(((int16_t*)&wrs[i]) + 1, sizeof(FIXP_WTP), vl);

        vtmp0 = vwmul_vv_i32m2(vtimelo, vre, vl);
        vtmp1 = vwmul_vv_i32m2(vtimehi, vim, vl);

        vmdct = vadd_vv_i32m2(vtmp0, vtmp1, vl);
        vmdct = vneg_v_i32m2(vmdct, vl);

        vsxw_v_i32m2(mdctData, vmdctidx, vmdct, vl);

        vtimeidx = vsub_vx_u16m1(vtimeidx, vl * 2, vl);
        vmdctidx = vsub_vx_u32m2(vmdctidx, vl * 4, vl);
      }
    }

    /* We pass the shortened folded data (-D-Cr,A-Br) to the MDCT function */
    dct_IV(mdctData, tl, &mdctData_e);

    pMdctData_e[n] = (SHORT)mdctData_e;

    timeData += tl;
    mdctData += tl;

    hMdct->prev_wrs = wrs;
    hMdct->prev_fr = fr;
    hMdct->prev_tl = tl;
  }

  return nSpec * tl;
}
#endif /* FUNCTION_mdct_block */

#endif /* defined(__riscv_xthead) && __riscv_v == 7000 */