/* Include files */

#include "modelInterface.h"
#include "m_UIzktXqQH1UJHNoLKJcwx.h"
#include <string.h>
#include "mwmathutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 208,   /* lineNo */
  "CarrierSynchronizer",               /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/CarrierSynchronizer.m"/* pathName */
};

static emlrtRSInfo b_emlrtRSI = { 1,   /* lineNo */
  "Helper",                            /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/compiled/+comm/+internal/Helper.p"/* pathName */
};

static emlrtRSInfo c_emlrtRSI = { 1,   /* lineNo */
  "System",                            /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/shared/system/coder/+matlab/+system/+coder/System.p"/* pathName */
};

static emlrtRSInfo d_emlrtRSI = { 1,   /* lineNo */
  "SystemProp",                        /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/shared/system/coder/+matlab/+system/+coder/SystemProp.p"/* pathName */
};

static emlrtRSInfo e_emlrtRSI = { 1,   /* lineNo */
  "SystemCore",                        /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/shared/system/coder/+matlab/+system/+coder/SystemCore.p"/* pathName */
};

static emlrtRSInfo f_emlrtRSI = { 12,  /* lineNo */
  "",                                  /* fcnName */
  ""                                   /* pathName */
};

static emlrtRSInfo g_emlrtRSI = { 22,  /* lineNo */
  "",                                  /* fcnName */
  ""                                   /* pathName */
};

static emlrtRSInfo h_emlrtRSI = { 23,  /* lineNo */
  "",                                  /* fcnName */
  ""                                   /* pathName */
};

static emlrtRSInfo i_emlrtRSI = { 24,  /* lineNo */
  "",                                  /* fcnName */
  ""                                   /* pathName */
};

static emlrtRSInfo j_emlrtRSI = { 25,  /* lineNo */
  "",                                  /* fcnName */
  ""                                   /* pathName */
};

static emlrtRSInfo k_emlrtRSI = { 31,  /* lineNo */
  "",                                  /* fcnName */
  ""                                   /* pathName */
};

static emlrtRSInfo l_emlrtRSI = { 236, /* lineNo */
  "CarrierSynchronizer",               /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/CarrierSynchronizer.m"/* pathName */
};

static emlrtRSInfo m_emlrtRSI = { 93,  /* lineNo */
  "validateattributes",                /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/lib/matlab/lang/validateattributes.m"/* pathName */
};

static emlrtRSInfo n_emlrtRSI = { 215, /* lineNo */
  "CarrierSynchronizer",               /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/CarrierSynchronizer.m"/* pathName */
};

static emlrtRSInfo o_emlrtRSI = { 222, /* lineNo */
  "CarrierSynchronizer",               /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/CarrierSynchronizer.m"/* pathName */
};

static emlrtRSInfo p_emlrtRSI = { 302, /* lineNo */
  "CarrierSynchronizer",               /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/CarrierSynchronizer.m"/* pathName */
};

static emlrtRSInfo q_emlrtRSI = { 33,  /* lineNo */
  "",                                  /* fcnName */
  ""                                   /* pathName */
};

static emlrtRSInfo r_emlrtRSI = { 39,  /* lineNo */
  "",                                  /* fcnName */
  ""                                   /* pathName */
};

static emlrtRSInfo s_emlrtRSI = { 48,  /* lineNo */
  "",                                  /* fcnName */
  ""                                   /* pathName */
};

static emlrtRSInfo t_emlrtRSI = { 57,  /* lineNo */
  "",                                  /* fcnName */
  ""                                   /* pathName */
};

static emlrtRSInfo u_emlrtRSI = { 65,  /* lineNo */
  "",                                  /* fcnName */
  ""                                   /* pathName */
};

static emlrtRSInfo v_emlrtRSI = { 412, /* lineNo */
  "CarrierSynchronizer",               /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/CarrierSynchronizer.m"/* pathName */
};

static emlrtRSInfo w_emlrtRSI = { 334, /* lineNo */
  "CarrierSynchronizer",               /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/CarrierSynchronizer.m"/* pathName */
};

static emlrtRSInfo x_emlrtRSI = { 34,  /* lineNo */
  "carrierSyncCore",                   /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/+internal/+carrierSync/carrierSyncCore.m"/* pathName */
};

static emlrtRSInfo y_emlrtRSI = { 35,  /* lineNo */
  "carrierSyncCore",                   /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/+internal/+carrierSync/carrierSyncCore.m"/* pathName */
};

static emlrtRSInfo ab_emlrtRSI = { 38, /* lineNo */
  "carrierSyncCore",                   /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/+internal/+carrierSync/carrierSyncCore.m"/* pathName */
};

static emlrtRSInfo bb_emlrtRSI = { 138,/* lineNo */
  "carrierSyncCore",                   /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/+internal/+carrierSync/carrierSyncCore.m"/* pathName */
};

static emlrtRSInfo cb_emlrtRSI = { 167,/* lineNo */
  "carrierSyncCore",                   /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/+internal/+carrierSync/carrierSyncCore.m"/* pathName */
};

static emlrtMCInfo emlrtMCI = { 14,    /* lineNo */
  37,                                  /* colNo */
  "validatepositive",                  /* fName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+valattr/validatepositive.m"/* pName */
};

static emlrtMCInfo b_emlrtMCI = { 14,  /* lineNo */
  37,                                  /* colNo */
  "validatefinite",                    /* fName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+valattr/validatefinite.m"/* pName */
};

static emlrtMCInfo c_emlrtMCI = { 14,  /* lineNo */
  37,                                  /* colNo */
  "validatenonnan",                    /* fName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+valattr/validatenonnan.m"/* pName */
};

static emlrtMCInfo d_emlrtMCI = { 13,  /* lineNo */
  37,                                  /* colNo */
  "validateinteger",                   /* fName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+valattr/validateinteger.m"/* pName */
};

static emlrtMCInfo e_emlrtMCI = { 22,  /* lineNo */
  27,                                  /* colNo */
  "validategt",                        /* fName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+valattr/validategt.m"/* pName */
};

static emlrtMCInfo f_emlrtMCI = { 22,  /* lineNo */
  27,                                  /* colNo */
  "validatele",                        /* fName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+valattr/validatele.m"/* pName */
};

static emlrtMCInfo g_emlrtMCI = { 1,   /* lineNo */
  1,                                   /* colNo */
  "SystemCore",                        /* fName */
  "/Applications/MATLAB_R2024a.app/toolbox/shared/system/coder/+matlab/+system/+coder/SystemCore.p"/* pName */
};

static emlrtECInfo emlrtECI = { 1,     /* nDims */
  3,                                   /* lineNo */
  4,                                   /* colNo */
  "",                                  /* fName */
  ""                                   /* pName */
};

static emlrtBCInfo emlrtBCI = { -1,    /* iFirst */
  -1,                                  /* iLast */
  38,                                  /* lineNo */
  20,                                  /* colNo */
  "",                                  /* aName */
  "carrierSyncCore",                   /* fName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/+internal/+carrierSync/carrierSyncCore.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo b_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  51,                                  /* lineNo */
  29,                                  /* colNo */
  "",                                  /* aName */
  "carrierSyncCore",                   /* fName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/+internal/+carrierSync/carrierSyncCore.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo c_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  52,                                  /* lineNo */
  37,                                  /* colNo */
  "",                                  /* aName */
  "carrierSyncCore",                   /* fName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/+internal/+carrierSync/carrierSyncCore.m",/* pName */
  0                                    /* checkKind */
};

static emlrtBCInfo d_emlrtBCI = { -1,  /* iFirst */
  -1,                                  /* iLast */
  38,                                  /* lineNo */
  25,                                  /* colNo */
  "",                                  /* aName */
  "carrierSyncCore",                   /* fName */
  "/Applications/MATLAB_R2024a.app/toolbox/comm/comm/+comm/+internal/+carrierSync/carrierSyncCore.m",/* pName */
  0                                    /* checkKind */
};

static emlrtRSInfo db_emlrtRSI = { 3,  /* lineNo */
  "",                                  /* fcnName */
  ""                                   /* pathName */
};

static emlrtRSInfo eb_emlrtRSI = { 14, /* lineNo */
  "validatenonnan",                    /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+valattr/validatenonnan.m"/* pathName */
};

static emlrtRSInfo fb_emlrtRSI = { 14, /* lineNo */
  "validatefinite",                    /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+valattr/validatefinite.m"/* pathName */
};

static emlrtRSInfo gb_emlrtRSI = { 14, /* lineNo */
  "validatepositive",                  /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+valattr/validatepositive.m"/* pathName */
};

static emlrtRSInfo hb_emlrtRSI = { 13, /* lineNo */
  "validateinteger",                   /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+valattr/validateinteger.m"/* pathName */
};

static emlrtRSInfo ib_emlrtRSI = { 22, /* lineNo */
  "validatele",                        /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+valattr/validatele.m"/* pathName */
};

static emlrtRSInfo jb_emlrtRSI = { 22, /* lineNo */
  "validategt",                        /* fcnName */
  "/Applications/MATLAB_R2024a.app/toolbox/eml/eml/+coder/+internal/+valattr/validategt.m"/* pathName */
};

/* Function Declarations */
static void cgxe_mdl_start(InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance);
static void cgxe_mdl_initialize(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance);
static void cgxe_mdl_outputs(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance);
static void cgxe_mdl_update(InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance);
static void cgxe_mdl_derivative(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance);
static void cgxe_mdl_enable(InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance);
static void cgxe_mdl_disable(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance);
static void cgxe_mdl_terminate(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance);
static void mw__internal__call__setup(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const emlrtStack *sp, real_T SamplesPerSymbol, real_T
  DampingFactor, real_T NormalizedLoopBandwidth);
static comm_CarrierSynchronizer *CarrierSynchronizer_CarrierSynchronizer
  (comm_CarrierSynchronizer *obj);
static void CarrierSynchronizer_set_SamplesPerSymbol(const emlrtStack *sp,
  comm_CarrierSynchronizer *obj, real_T b_value);
static void CarrierSynchronizer_set_DampingFactor(const emlrtStack *sp,
  comm_CarrierSynchronizer *obj, real_T b_value);
static void CarrierSynchronizer_set_NormalizedLoopBandwidth(const emlrtStack *sp,
  comm_CarrierSynchronizer *obj, real_T b_value);
static void CarrierSynchronizer_CalculateLoopGains(comm_CarrierSynchronizer *obj);
static void SystemCore_checkTunablePropChange(const emlrtStack *sp,
  comm_CarrierSynchronizer *obj);
static void mw__internal__call__reset(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const emlrtStack *sp, real_T SamplesPerSymbol, real_T
  DampingFactor, real_T NormalizedLoopBandwidth);
static void mw__internal__call__step(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const emlrtStack *sp, real_T SamplesPerSymbol, real_T
  DampingFactor, real_T NormalizedLoopBandwidth, coder_array_creal_T_2D *u0,
  coder_array_creal_T_2D *b_y0);
static void mw__internal__system___fcn(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const emlrtStack *sp, real_T varargin_3, real_T varargin_4,
  real_T varargin_5, coder_array_creal_T *varargin_7, coder_array_creal_T
  *varargout_1);
static void CarrierSynchronizer_stepImpl(const emlrtStack *sp,
  comm_CarrierSynchronizer *obj, coder_array_creal_T *input, coder_array_creal_T
  *output);
static creal_T expFunc(real_T in);
static const mxArray *emlrt_marshallOut(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const emlrtStack *sp);
static const mxArray *cgxe_mdl_get_sim_state
  (InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance);
static void emlrt_marshallIn(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const emlrtStack *sp, const mxArray *u);
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
  const char_T *identifier, int32_T y[2]);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, int32_T y[2]);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
  const char_T *identifier, comm_CarrierSynchronizer *y);
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, comm_CarrierSynchronizer *y);
static int32_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static boolean_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, cell_wrap y[1]);
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, uint32_T y[8]);
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[4]);
static real_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static creal_T l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId);
static boolean_T m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
  const char_T *identifier);
static void cgxe_mdl_set_sim_state(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const mxArray *st);
static const mxArray *message(const emlrtStack *sp, const mxArray *m1, const
  mxArray *m2, emlrtMCInfo *location);
static const mxArray *getString(const emlrtStack *sp, const mxArray *m1,
  emlrtMCInfo *location);
static void error(const emlrtStack *sp, const mxArray *m, const mxArray *m1,
                  emlrtMCInfo *location);
static const mxArray *b_message(const emlrtStack *sp, const mxArray *m1, const
  mxArray *m2, const mxArray *m3, const mxArray *m4, emlrtMCInfo *location);
static const mxArray *c_message(const emlrtStack *sp, const mxArray *m1,
  emlrtMCInfo *location);
static void n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, int32_T ret[2]);
static int32_T o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId);
static boolean_T p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId);
static void q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, uint32_T ret[8]);
static void r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[4]);
static real_T s_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId);
static creal_T t_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId);
static void carrierSyncCore(const emlrtStack *sp, coder_array_creal_T *input,
  real_T *loopFiltState, real_T *integFiltState, real_T *DDSPreviousInp, creal_T
  *previousSample, real_T *phase, real_T *digitalSynthesizerGain, real_T
  *integratorGain, real_T *proportionalGain, real_T *phaseOffset,
  coder_array_creal_T *output, real_T phaseEstimate_data[], int32_T
  phaseEstimate_size[1]);
static void array_creal_T_2D_SetSize(coder_array_creal_T_2D *coderArray, int32_T
  size0, int32_T size1);
static void array_creal_T_SetSize(coder_array_creal_T *coderArray, int32_T size0);
static void array_creal_T_2D_Constructor(coder_array_creal_T_2D *coderArray);
static void array_creal_T_2D_Destructor(coder_array_creal_T_2D *coderArray);
static void array_creal_T_Constructor(coder_array_creal_T *coderArray);
static void array_creal_T_Destructor(coder_array_creal_T *coderArray);
static void init_simulink_io_address(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance);

/* Function Definitions */
static void cgxe_mdl_start(InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  real_T *DampingFactor;
  real_T *NormalizedLoopBandwidth;
  real_T *SamplesPerSymbol;
  init_simulink_io_address(moduleInstance);
  NormalizedLoopBandwidth = (real_T *)cgxertGetRunTimeParamInfoData
    (moduleInstance->S, 2);
  DampingFactor = (real_T *)cgxertGetRunTimeParamInfoData(moduleInstance->S, 1);
  SamplesPerSymbol = (real_T *)cgxertGetRunTimeParamInfoData(moduleInstance->S,
    0);
  st.tls = moduleInstance->emlrtRootTLSGlobal;
  cgxertSetGcb(moduleInstance->S, -1, -1);
  mw__internal__call__setup(moduleInstance, &st, *SamplesPerSymbol,
    *DampingFactor, *NormalizedLoopBandwidth);
  cgxertRestoreGcb(moduleInstance->S, -1, -1);
}

static void cgxe_mdl_initialize(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance)
{
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  real_T *DampingFactor;
  real_T *NormalizedLoopBandwidth;
  real_T *SamplesPerSymbol;
  NormalizedLoopBandwidth = (real_T *)cgxertGetRunTimeParamInfoData
    (moduleInstance->S, 2);
  DampingFactor = (real_T *)cgxertGetRunTimeParamInfoData(moduleInstance->S, 1);
  SamplesPerSymbol = (real_T *)cgxertGetRunTimeParamInfoData(moduleInstance->S,
    0);
  st.tls = moduleInstance->emlrtRootTLSGlobal;
  emlrtLicenseCheckR2022a(&st, "EMLRT:runTime:MexFunctionNeedsLicense",
    "communication_toolbox", 2);
  cgxertSetGcb(moduleInstance->S, -1, -1);
  mw__internal__call__reset(moduleInstance, &st, *SamplesPerSymbol,
    *DampingFactor, *NormalizedLoopBandwidth);
  cgxertRestoreGcb(moduleInstance->S, -1, -1);
}

static void cgxe_mdl_outputs(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance)
{
  coder_array_creal_T_2D u_tmp3;
  coder_array_creal_T_2D y_tmp0;
  emlrtStack st = { NULL,              /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  real_T *DampingFactor;
  real_T *NormalizedLoopBandwidth;
  real_T *SamplesPerSymbol;
  int32_T i;
  int32_T loop_ub;
  NormalizedLoopBandwidth = (real_T *)cgxertGetRunTimeParamInfoData
    (moduleInstance->S, 2);
  DampingFactor = (real_T *)cgxertGetRunTimeParamInfoData(moduleInstance->S, 1);
  SamplesPerSymbol = (real_T *)cgxertGetRunTimeParamInfoData(moduleInstance->S,
    0);
  st.tls = moduleInstance->emlrtRootTLSGlobal;
  cgxertSetGcb(moduleInstance->S, -1, -1);
  cgxertCheckCurrentInputPortDimensions(moduleInstance->S,
    "comm.CarrierSynchronizer", "stepImpl");
  array_creal_T_2D_Constructor(&u_tmp3);
  array_creal_T_2D_SetSize(&u_tmp3, (*moduleInstance->u0_sizes)[0],
    (*moduleInstance->u0_sizes)[1]);
  loop_ub = (*moduleInstance->u0_sizes)[0] * (*moduleInstance->u0_sizes)[1];
  for (i = 0; i < loop_ub; i++) {
    u_tmp3.vector.data[i].re = (*moduleInstance->u0_data)[i].re;
    u_tmp3.vector.data[i].im = (*moduleInstance->u0_data)[i].im;
  }

  array_creal_T_2D_Constructor(&y_tmp0);
  mw__internal__call__step(moduleInstance, &st, *SamplesPerSymbol,
    *DampingFactor, *NormalizedLoopBandwidth, &u_tmp3, &y_tmp0);
  array_creal_T_2D_Destructor(&u_tmp3);
  cgxertSetCurrentOutputPortDimensions(moduleInstance->S, 0, 0, y_tmp0.size[0]);
  cgxertSetCurrentOutputPortDimensions(moduleInstance->S, 0, 1, 1);
  if (y_tmp0.size[0] - 1 >= 0) {
    memcpy(&(*moduleInstance->y0_data)[0], &y_tmp0.vector.data[0], (uint32_T)
           y_tmp0.size[0] * sizeof(creal_T));
  }

  array_creal_T_2D_Destructor(&y_tmp0);
  cgxertRestoreGcb(moduleInstance->S, -1, -1);
}

static void cgxe_mdl_update(InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance)
{
  (void)moduleInstance;
}

static void cgxe_mdl_derivative(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance)
{
  (void)moduleInstance;
}

static void cgxe_mdl_enable(InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance)
{
  (void)moduleInstance;
}

static void cgxe_mdl_disable(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance)
{
  (void)moduleInstance;
}

static void cgxe_mdl_terminate(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance)
{
  cgxertSetGcb(moduleInstance->S, -1, -1);
  cgxertRestoreGcb(moduleInstance->S, -1, -1);
}

static void mw__internal__call__setup(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const emlrtStack *sp, real_T SamplesPerSymbol, real_T
  DampingFactor, real_T NormalizedLoopBandwidth)
{
  static const int32_T iv[2] = { 1, 51 };

  static const int32_T iv2[2] = { 1, 51 };

  static const int32_T iv3[2] = { 1, 5 };

  static int16_T iv1[8] = { 6175, 1, 1, 1, 1, 1, 1, 1 };

  static char_T b_u[51] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 's', 'y', 's',
    't', 'e', 'm', ':', 'm', 'e', 't', 'h', 'o', 'd', 'C', 'a', 'l', 'l', 'e',
    'd', 'W', 'h', 'e', 'n', 'L', 'o', 'c', 'k', 'e', 'd', 'R', 'e', 'l', 'e',
    'a', 's', 'e', 'd', 'C', 'o', 'd', 'e', 'g', 'e', 'n' };

  static char_T d_u[5] = { 's', 'e', 't', 'u', 'p' };

  static char_T correctedValue[4] = { 'A', 'u', 't', 'o' };

  cell_wrap varSizes[1];
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *m;
  const mxArray *y;
  int32_T i;
  char_T u[51];
  char_T c_u[5];
  boolean_T flag;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  if (!moduleInstance->sysobj_not_empty) {
    st.site = &f_emlrtRSI;
    CarrierSynchronizer_CarrierSynchronizer(&moduleInstance->sysobj);
    moduleInstance->sysobj_not_empty = true;
    st.site = &g_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    c_st.site = &d_emlrtRSI;
    flag = (moduleInstance->sysobj.isInitialized == 1);
    if (flag) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    for (i = 0; i < 4; i++) {
      moduleInstance->sysobj.ModulationPhaseOffset[i] = correctedValue[i];
    }

    st.site = &h_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    flag = (moduleInstance->sysobj.isInitialized == 1);
    if (flag) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    st.site = &h_emlrtRSI;
    CarrierSynchronizer_set_SamplesPerSymbol(&st, &moduleInstance->sysobj,
      SamplesPerSymbol);
    st.site = &i_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    flag = (moduleInstance->sysobj.isInitialized == 1);
    if (flag) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    st.site = &i_emlrtRSI;
    CarrierSynchronizer_set_DampingFactor(&st, &moduleInstance->sysobj,
      DampingFactor);
    st.site = &j_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    flag = (moduleInstance->sysobj.isInitialized == 1);
    if (flag) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    st.site = &j_emlrtRSI;
    CarrierSynchronizer_set_NormalizedLoopBandwidth(&st, &moduleInstance->sysobj,
      NormalizedLoopBandwidth);
  }

  st.site = &k_emlrtRSI;
  if (moduleInstance->sysobj.isInitialized != 0) {
    for (i = 0; i < 51; i++) {
      u[i] = b_u[i];
    }

    y = NULL;
    m = emlrtCreateCharArray(2, &iv[0]);
    emlrtInitCharArrayR2013a(&st, 51, m, &u[0]);
    emlrtAssign(&y, m);
    for (i = 0; i < 51; i++) {
      u[i] = b_u[i];
    }

    b_y = NULL;
    m = emlrtCreateCharArray(2, &iv2[0]);
    emlrtInitCharArrayR2013a(&st, 51, m, &u[0]);
    emlrtAssign(&b_y, m);
    for (i = 0; i < 5; i++) {
      c_u[i] = d_u[i];
    }

    c_y = NULL;
    m = emlrtCreateCharArray(2, &iv3[0]);
    emlrtInitCharArrayR2013a(&st, 5, m, &c_u[0]);
    emlrtAssign(&c_y, m);
    b_st.site = &e_emlrtRSI;
    error(&b_st, y, getString(&b_st, message(&b_st, b_y, c_y, &g_emlrtMCI),
           &g_emlrtMCI), &g_emlrtMCI);
  }

  moduleInstance->sysobj.isInitialized = 1;
  b_st.site = &e_emlrtRSI;
  for (i = 0; i < 8; i++) {
    varSizes[0].f1[i] = (uint32_T)iv1[i];
  }

  moduleInstance->sysobj.inputVarSize[0] = varSizes[0];
  b_st.site = &e_emlrtRSI;
  moduleInstance->sysobj.pPhase = 0.0;
  moduleInstance->sysobj.pPreviousSample.re = 0.0;
  moduleInstance->sysobj.pPreviousSample.im = 0.0;
  c_st.site = &p_emlrtRSI;
  CarrierSynchronizer_CalculateLoopGains(&moduleInstance->sysobj);
  moduleInstance->sysobj.pDigitalSynthesizerGain = -1.0;
  b_st.site = &e_emlrtRSI;
  SystemCore_checkTunablePropChange(&b_st, &moduleInstance->sysobj);
  moduleInstance->sysobj.TunablePropsChanged = false;
}

static comm_CarrierSynchronizer *CarrierSynchronizer_CarrierSynchronizer
  (comm_CarrierSynchronizer *obj)
{
  comm_CarrierSynchronizer *b_obj;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack st;
  b_obj = obj;
  b_obj->CustomPhaseOffset = 0.0;
  st.site = &emlrtRSI;
  b_st.site = &b_emlrtRSI;
  c_st.site = &c_emlrtRSI;
  d_st.site = &d_emlrtRSI;
  c_st.site = &c_emlrtRSI;
  b_obj->TunablePropsChanged = false;
  b_obj->CacheInputSizes = false;
  d_st.site = &e_emlrtRSI;
  b_obj->isInitialized = 0;
  return b_obj;
}

static void CarrierSynchronizer_set_SamplesPerSymbol(const emlrtStack *sp,
  comm_CarrierSynchronizer *obj, real_T b_value)
{
  static const int32_T iv[2] = { 1, 23 };

  static const int32_T iv1[2] = { 1, 21 };

  static const int32_T iv10[2] = { 1, 16 };

  static const int32_T iv11[2] = { 1, 16 };

  static const int32_T iv2[2] = { 1, 48 };

  static const int32_T iv3[2] = { 1, 21 };

  static const int32_T iv4[2] = { 1, 46 };

  static const int32_T iv5[2] = { 1, 16 };

  static const int32_T iv6[2] = { 1, 22 };

  static const int32_T iv7[2] = { 1, 46 };

  static const int32_T iv8[2] = { 1, 16 };

  static const int32_T iv9[2] = { 1, 47 };

  static char_T f_u[48] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'V', 'a', 'l', 'i', 'd', 'a', 't', 'e', 'a', 't', 't',
    'r', 'i', 'b', 'u', 't', 'e', 's', 'e', 'x', 'p', 'e', 'c', 't', 'e', 'd',
    'P', 'o', 's', 'i', 't', 'i', 'v', 'e' };

  static char_T p_u[47] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'V', 'a', 'l', 'i', 'd', 'a', 't', 'e', 'a', 't', 't',
    'r', 'i', 'b', 'u', 't', 'e', 's', 'e', 'x', 'p', 'e', 'c', 't', 'e', 'd',
    'I', 'n', 't', 'e', 'g', 'e', 'r' };

  static char_T j_u[46] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'V', 'a', 'l', 'i', 'd', 'a', 't', 'e', 'a', 't', 't',
    'r', 'i', 'b', 'u', 't', 'e', 's', 'e', 'x', 'p', 'e', 'c', 't', 'e', 'd',
    'F', 'i', 'n', 'i', 't', 'e' };

  static char_T n_u[46] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'V', 'a', 'l', 'i', 'd', 'a', 't', 'e', 'a', 't', 't',
    'r', 'i', 'b', 'u', 't', 'e', 's', 'e', 'x', 'p', 'e', 'c', 't', 'e', 'd',
    'N', 'o', 'n', 'N', 'a', 'N' };

  static char_T b_u[23] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'e', 'x', 'p',
    'e', 'c', 't', 'e', 'd', 'P', 'o', 's', 'i', 't', 'i', 'v', 'e' };

  static char_T m_u[22] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'e', 'x', 'p',
    'e', 'c', 't', 'e', 'd', 'I', 'n', 't', 'e', 'g', 'e', 'r' };

  static char_T e_u[21] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'e', 'x', 'p',
    'e', 'c', 't', 'e', 'd', 'F', 'i', 'n', 'i', 't', 'e' };

  static char_T i_u[21] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'e', 'x', 'p',
    'e', 'c', 't', 'e', 'd', 'N', 'o', 'n', 'N', 'a', 'N' };

  static char_T k_u[16] = { 'S', 'a', 'm', 'p', 'l', 'e', 's', 'P', 'e', 'r',
    'S', 'y', 'm', 'b', 'o', 'l' };

  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *f_y;
  const mxArray *g_y;
  const mxArray *h_y;
  const mxArray *i_y;
  const mxArray *j_y;
  const mxArray *k_y;
  const mxArray *l_y;
  const mxArray *m;
  const mxArray *y;
  int32_T i;
  char_T d_u[48];
  char_T o_u[47];
  char_T g_u[46];
  char_T u[23];
  char_T l_u[22];
  char_T c_u[21];
  char_T h_u[16];
  boolean_T p;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &l_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  b_st.site = &m_emlrtRSI;
  p = true;
  if (b_value <= 0.0) {
    p = false;
  }

  if (!p) {
    for (i = 0; i < 23; i++) {
      u[i] = b_u[i];
    }

    y = NULL;
    m = emlrtCreateCharArray(2, &iv[0]);
    emlrtInitCharArrayR2013a(&b_st, 23, m, &u[0]);
    emlrtAssign(&y, m);
    for (i = 0; i < 48; i++) {
      d_u[i] = f_u[i];
    }

    c_y = NULL;
    m = emlrtCreateCharArray(2, &iv2[0]);
    emlrtInitCharArrayR2013a(&b_st, 48, m, &d_u[0]);
    emlrtAssign(&c_y, m);
    for (i = 0; i < 16; i++) {
      h_u[i] = k_u[i];
    }

    f_y = NULL;
    m = emlrtCreateCharArray(2, &iv5[0]);
    emlrtInitCharArrayR2013a(&b_st, 16, m, &h_u[0]);
    emlrtAssign(&f_y, m);
    c_st.site = &gb_emlrtRSI;
    error(&c_st, y, getString(&c_st, message(&c_st, c_y, f_y, &emlrtMCI),
           &emlrtMCI), &emlrtMCI);
  }

  b_st.site = &m_emlrtRSI;
  p = true;
  if ((!!muDoubleScalarIsInf(b_value)) || (!!muDoubleScalarIsNaN(b_value))) {
    p = false;
  }

  if (!p) {
    for (i = 0; i < 21; i++) {
      c_u[i] = e_u[i];
    }

    b_y = NULL;
    m = emlrtCreateCharArray(2, &iv1[0]);
    emlrtInitCharArrayR2013a(&b_st, 21, m, &c_u[0]);
    emlrtAssign(&b_y, m);
    for (i = 0; i < 46; i++) {
      g_u[i] = j_u[i];
    }

    e_y = NULL;
    m = emlrtCreateCharArray(2, &iv4[0]);
    emlrtInitCharArrayR2013a(&b_st, 46, m, &g_u[0]);
    emlrtAssign(&e_y, m);
    for (i = 0; i < 16; i++) {
      h_u[i] = k_u[i];
    }

    i_y = NULL;
    m = emlrtCreateCharArray(2, &iv8[0]);
    emlrtInitCharArrayR2013a(&b_st, 16, m, &h_u[0]);
    emlrtAssign(&i_y, m);
    c_st.site = &fb_emlrtRSI;
    error(&c_st, b_y, getString(&c_st, message(&c_st, e_y, i_y, &b_emlrtMCI),
           &b_emlrtMCI), &b_emlrtMCI);
  }

  b_st.site = &m_emlrtRSI;
  p = true;
  if (muDoubleScalarIsNaN(b_value)) {
    p = false;
  }

  if (!p) {
    for (i = 0; i < 21; i++) {
      c_u[i] = i_u[i];
    }

    d_y = NULL;
    m = emlrtCreateCharArray(2, &iv3[0]);
    emlrtInitCharArrayR2013a(&b_st, 21, m, &c_u[0]);
    emlrtAssign(&d_y, m);
    for (i = 0; i < 46; i++) {
      g_u[i] = n_u[i];
    }

    h_y = NULL;
    m = emlrtCreateCharArray(2, &iv7[0]);
    emlrtInitCharArrayR2013a(&b_st, 46, m, &g_u[0]);
    emlrtAssign(&h_y, m);
    for (i = 0; i < 16; i++) {
      h_u[i] = k_u[i];
    }

    k_y = NULL;
    m = emlrtCreateCharArray(2, &iv10[0]);
    emlrtInitCharArrayR2013a(&b_st, 16, m, &h_u[0]);
    emlrtAssign(&k_y, m);
    c_st.site = &eb_emlrtRSI;
    error(&c_st, d_y, getString(&c_st, message(&c_st, h_y, k_y, &c_emlrtMCI),
           &c_emlrtMCI), &c_emlrtMCI);
  }

  b_st.site = &m_emlrtRSI;
  p = true;
  if ((!muDoubleScalarIsInf(b_value)) && (!muDoubleScalarIsNaN(b_value)) &&
      (muDoubleScalarFloor(b_value) == b_value)) {
  } else {
    p = false;
  }

  if (!p) {
    for (i = 0; i < 22; i++) {
      l_u[i] = m_u[i];
    }

    g_y = NULL;
    m = emlrtCreateCharArray(2, &iv6[0]);
    emlrtInitCharArrayR2013a(&b_st, 22, m, &l_u[0]);
    emlrtAssign(&g_y, m);
    for (i = 0; i < 47; i++) {
      o_u[i] = p_u[i];
    }

    j_y = NULL;
    m = emlrtCreateCharArray(2, &iv9[0]);
    emlrtInitCharArrayR2013a(&b_st, 47, m, &o_u[0]);
    emlrtAssign(&j_y, m);
    for (i = 0; i < 16; i++) {
      h_u[i] = k_u[i];
    }

    l_y = NULL;
    m = emlrtCreateCharArray(2, &iv11[0]);
    emlrtInitCharArrayR2013a(&b_st, 16, m, &h_u[0]);
    emlrtAssign(&l_y, m);
    c_st.site = &hb_emlrtRSI;
    error(&c_st, g_y, getString(&c_st, message(&c_st, j_y, l_y, &d_emlrtMCI),
           &d_emlrtMCI), &d_emlrtMCI);
  }

  obj->SamplesPerSymbol = b_value;
}

static void CarrierSynchronizer_set_DampingFactor(const emlrtStack *sp,
  comm_CarrierSynchronizer *obj, real_T b_value)
{
  static const int32_T iv[2] = { 1, 23 };

  static const int32_T iv1[2] = { 1, 21 };

  static const int32_T iv2[2] = { 1, 48 };

  static const int32_T iv3[2] = { 1, 21 };

  static const int32_T iv4[2] = { 1, 46 };

  static const int32_T iv5[2] = { 1, 13 };

  static const int32_T iv6[2] = { 1, 46 };

  static const int32_T iv7[2] = { 1, 13 };

  static const int32_T iv8[2] = { 1, 13 };

  static char_T f_u[48] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'V', 'a', 'l', 'i', 'd', 'a', 't', 'e', 'a', 't', 't',
    'r', 'i', 'b', 'u', 't', 'e', 's', 'e', 'x', 'p', 'e', 'c', 't', 'e', 'd',
    'P', 'o', 's', 'i', 't', 'i', 'v', 'e' };

  static char_T j_u[46] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'V', 'a', 'l', 'i', 'd', 'a', 't', 'e', 'a', 't', 't',
    'r', 'i', 'b', 'u', 't', 'e', 's', 'e', 'x', 'p', 'e', 'c', 't', 'e', 'd',
    'F', 'i', 'n', 'i', 't', 'e' };

  static char_T l_u[46] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'V', 'a', 'l', 'i', 'd', 'a', 't', 'e', 'a', 't', 't',
    'r', 'i', 'b', 'u', 't', 'e', 's', 'e', 'x', 'p', 'e', 'c', 't', 'e', 'd',
    'N', 'o', 'n', 'N', 'a', 'N' };

  static char_T b_u[23] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'e', 'x', 'p',
    'e', 'c', 't', 'e', 'd', 'P', 'o', 's', 'i', 't', 'i', 'v', 'e' };

  static char_T e_u[21] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'e', 'x', 'p',
    'e', 'c', 't', 'e', 'd', 'F', 'i', 'n', 'i', 't', 'e' };

  static char_T i_u[21] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'e', 'x', 'p',
    'e', 'c', 't', 'e', 'd', 'N', 'o', 'n', 'N', 'a', 'N' };

  static char_T k_u[13] = { 'D', 'a', 'm', 'p', 'i', 'n', 'g', 'F', 'a', 'c',
    't', 'o', 'r' };

  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *f_y;
  const mxArray *g_y;
  const mxArray *h_y;
  const mxArray *i_y;
  const mxArray *m;
  const mxArray *y;
  int32_T i;
  char_T d_u[48];
  char_T g_u[46];
  char_T u[23];
  char_T c_u[21];
  char_T h_u[13];
  boolean_T p;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &n_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  b_st.site = &m_emlrtRSI;
  p = true;
  if (b_value <= 0.0) {
    p = false;
  }

  if (!p) {
    for (i = 0; i < 23; i++) {
      u[i] = b_u[i];
    }

    y = NULL;
    m = emlrtCreateCharArray(2, &iv[0]);
    emlrtInitCharArrayR2013a(&b_st, 23, m, &u[0]);
    emlrtAssign(&y, m);
    for (i = 0; i < 48; i++) {
      d_u[i] = f_u[i];
    }

    c_y = NULL;
    m = emlrtCreateCharArray(2, &iv2[0]);
    emlrtInitCharArrayR2013a(&b_st, 48, m, &d_u[0]);
    emlrtAssign(&c_y, m);
    for (i = 0; i < 13; i++) {
      h_u[i] = k_u[i];
    }

    f_y = NULL;
    m = emlrtCreateCharArray(2, &iv5[0]);
    emlrtInitCharArrayR2013a(&b_st, 13, m, &h_u[0]);
    emlrtAssign(&f_y, m);
    c_st.site = &gb_emlrtRSI;
    error(&c_st, y, getString(&c_st, message(&c_st, c_y, f_y, &emlrtMCI),
           &emlrtMCI), &emlrtMCI);
  }

  b_st.site = &m_emlrtRSI;
  p = true;
  if ((!!muDoubleScalarIsInf(b_value)) || (!!muDoubleScalarIsNaN(b_value))) {
    p = false;
  }

  if (!p) {
    for (i = 0; i < 21; i++) {
      c_u[i] = e_u[i];
    }

    b_y = NULL;
    m = emlrtCreateCharArray(2, &iv1[0]);
    emlrtInitCharArrayR2013a(&b_st, 21, m, &c_u[0]);
    emlrtAssign(&b_y, m);
    for (i = 0; i < 46; i++) {
      g_u[i] = j_u[i];
    }

    e_y = NULL;
    m = emlrtCreateCharArray(2, &iv4[0]);
    emlrtInitCharArrayR2013a(&b_st, 46, m, &g_u[0]);
    emlrtAssign(&e_y, m);
    for (i = 0; i < 13; i++) {
      h_u[i] = k_u[i];
    }

    h_y = NULL;
    m = emlrtCreateCharArray(2, &iv7[0]);
    emlrtInitCharArrayR2013a(&b_st, 13, m, &h_u[0]);
    emlrtAssign(&h_y, m);
    c_st.site = &fb_emlrtRSI;
    error(&c_st, b_y, getString(&c_st, message(&c_st, e_y, h_y, &b_emlrtMCI),
           &b_emlrtMCI), &b_emlrtMCI);
  }

  b_st.site = &m_emlrtRSI;
  p = true;
  if (muDoubleScalarIsNaN(b_value)) {
    p = false;
  }

  if (!p) {
    for (i = 0; i < 21; i++) {
      c_u[i] = i_u[i];
    }

    d_y = NULL;
    m = emlrtCreateCharArray(2, &iv3[0]);
    emlrtInitCharArrayR2013a(&b_st, 21, m, &c_u[0]);
    emlrtAssign(&d_y, m);
    for (i = 0; i < 46; i++) {
      g_u[i] = l_u[i];
    }

    g_y = NULL;
    m = emlrtCreateCharArray(2, &iv6[0]);
    emlrtInitCharArrayR2013a(&b_st, 46, m, &g_u[0]);
    emlrtAssign(&g_y, m);
    for (i = 0; i < 13; i++) {
      h_u[i] = k_u[i];
    }

    i_y = NULL;
    m = emlrtCreateCharArray(2, &iv8[0]);
    emlrtInitCharArrayR2013a(&b_st, 13, m, &h_u[0]);
    emlrtAssign(&i_y, m);
    c_st.site = &eb_emlrtRSI;
    error(&c_st, d_y, getString(&c_st, message(&c_st, g_y, i_y, &c_emlrtMCI),
           &c_emlrtMCI), &c_emlrtMCI);
  }

  obj->DampingFactor = b_value;
}

static void CarrierSynchronizer_set_NormalizedLoopBandwidth(const emlrtStack *sp,
  comm_CarrierSynchronizer *obj, real_T b_value)
{
  static const int32_T iv[2] = { 1, 17 };

  static const int32_T iv1[2] = { 1, 19 };

  static const int32_T iv10[2] = { 1, 46 };

  static const int32_T iv11[2] = { 1, 23 };

  static const int32_T iv12[2] = { 1, 2 };

  static const int32_T iv13[2] = { 1, 46 };

  static const int32_T iv14[2] = { 1, 23 };

  static const int32_T iv15[2] = { 1, 23 };

  static const int32_T iv2[2] = { 1, 40 };

  static const int32_T iv3[2] = { 1, 23 };

  static const int32_T iv4[2] = { 1, 40 };

  static const int32_T iv5[2] = { 1, 23 };

  static const int32_T iv6[2] = { 1, 21 };

  static const int32_T iv7[2] = { 1, 48 };

  static const int32_T iv8[2] = { 1, 23 };

  static const int32_T iv9[2] = { 1, 21 };

  static char_T m_u[48] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'V', 'a', 'l', 'i', 'd', 'a', 't', 'e', 'a', 't', 't',
    'r', 'i', 'b', 'u', 't', 'e', 's', 'e', 'x', 'p', 'e', 'c', 't', 'e', 'd',
    'P', 'o', 's', 'i', 't', 'i', 'v', 'e' };

  static char_T q_u[46] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'V', 'a', 'l', 'i', 'd', 'a', 't', 'e', 'a', 't', 't',
    'r', 'i', 'b', 'u', 't', 'e', 's', 'e', 'x', 'p', 'e', 'c', 't', 'e', 'd',
    'F', 'i', 'n', 'i', 't', 'e' };

  static char_T s_u[46] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o', 'l',
    'b', 'o', 'x', ':', 'V', 'a', 'l', 'i', 'd', 'a', 't', 'e', 'a', 't', 't',
    'r', 'i', 'b', 'u', 't', 'e', 's', 'e', 'x', 'p', 'e', 'c', 't', 'e', 'd',
    'N', 'o', 'n', 'N', 'a', 'N' };

  static char_T f_u[40] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'v', 'a', 'l',
    'i', 'd', 'a', 't', 'e', 'a', 't', 't', 'r', 'i', 'b', 'u', 't', 'e', 's',
    ':', 'e', 'x', 'p', 'e', 'c', 't', 'e', 'd', 'S', 'c', 'a', 'l', 'a', 'r' };

  static char_T h_u[23] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'e', 'x', 'p',
    'e', 'c', 't', 'e', 'd', 'P', 'o', 's', 'i', 't', 'i', 'v', 'e' };

  static char_T i_u[23] = { 'N', 'o', 'r', 'm', 'a', 'l', 'i', 'z', 'e', 'd',
    'L', 'o', 'o', 'p', 'B', 'a', 'n', 'd', 'w', 'i', 'd', 't', 'h' };

  static char_T l_u[21] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'e', 'x', 'p',
    'e', 'c', 't', 'e', 'd', 'F', 'i', 'n', 'i', 't', 'e' };

  static char_T o_u[21] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'e', 'x', 'p',
    'e', 'c', 't', 'e', 'd', 'N', 'o', 'n', 'N', 'a', 'N' };

  static char_T e_u[19] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'n', 'o', 't',
    'L', 'e', 's', 's', 'E', 'q', 'u', 'a', 'l' };

  static char_T b_u[17] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'n', 'o', 't',
    'G', 'r', 'e', 'a', 't', 'e', 'r' };

  static char_T r_u[2] = { '<', '=' };

  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *f_y;
  const mxArray *g_y;
  const mxArray *h_y;
  const mxArray *i_y;
  const mxArray *j_y;
  const mxArray *k_y;
  const mxArray *l_y;
  const mxArray *m;
  const mxArray *m_y;
  const mxArray *n_y;
  const mxArray *o_y;
  const mxArray *p_y;
  const mxArray *q_y;
  const mxArray *r_y;
  const mxArray *s_y;
  const mxArray *y;
  int32_T i;
  char_T k_u[48];
  char_T n_u[46];
  char_T d_u[40];
  char_T g_u[23];
  char_T j_u[21];
  char_T c_u[19];
  char_T u[17];
  char_T p_u[2];
  boolean_T p;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &o_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  b_st.site = &m_emlrtRSI;
  p = true;
  if (!(b_value > 0.0)) {
    p = false;
  }

  if (!p) {
    for (i = 0; i < 17; i++) {
      u[i] = b_u[i];
    }

    y = NULL;
    m = emlrtCreateCharArray(2, &iv[0]);
    emlrtInitCharArrayR2013a(&b_st, 17, m, &u[0]);
    emlrtAssign(&y, m);
    for (i = 0; i < 40; i++) {
      d_u[i] = f_u[i];
    }

    c_y = NULL;
    m = emlrtCreateCharArray(2, &iv2[0]);
    emlrtInitCharArrayR2013a(&b_st, 40, m, &d_u[0]);
    emlrtAssign(&c_y, m);
    for (i = 0; i < 23; i++) {
      g_u[i] = i_u[i];
    }

    f_y = NULL;
    m = emlrtCreateCharArray(2, &iv5[0]);
    emlrtInitCharArrayR2013a(&b_st, 23, m, &g_u[0]);
    emlrtAssign(&f_y, m);
    h_y = NULL;
    m = emlrtCreateString1R2022a(&b_st, '>');
    emlrtAssign(&h_y, m);
    k_y = NULL;
    m = emlrtCreateString1R2022a(&b_st, '0');
    emlrtAssign(&k_y, m);
    c_st.site = &jb_emlrtRSI;
    error(&c_st, y, getString(&c_st, b_message(&c_st, c_y, f_y, h_y, k_y,
            &e_emlrtMCI), &e_emlrtMCI), &e_emlrtMCI);
  }

  b_st.site = &m_emlrtRSI;
  p = true;
  if (!(b_value <= 1.0)) {
    p = false;
  }

  if (!p) {
    for (i = 0; i < 19; i++) {
      c_u[i] = e_u[i];
    }

    b_y = NULL;
    m = emlrtCreateCharArray(2, &iv1[0]);
    emlrtInitCharArrayR2013a(&b_st, 19, m, &c_u[0]);
    emlrtAssign(&b_y, m);
    for (i = 0; i < 40; i++) {
      d_u[i] = f_u[i];
    }

    e_y = NULL;
    m = emlrtCreateCharArray(2, &iv4[0]);
    emlrtInitCharArrayR2013a(&b_st, 40, m, &d_u[0]);
    emlrtAssign(&e_y, m);
    for (i = 0; i < 23; i++) {
      g_u[i] = i_u[i];
    }

    j_y = NULL;
    m = emlrtCreateCharArray(2, &iv8[0]);
    emlrtInitCharArrayR2013a(&b_st, 23, m, &g_u[0]);
    emlrtAssign(&j_y, m);
    for (i = 0; i < 2; i++) {
      p_u[i] = r_u[i];
    }

    o_y = NULL;
    m = emlrtCreateCharArray(2, &iv12[0]);
    emlrtInitCharArrayR2013a(&b_st, 2, m, &p_u[0]);
    emlrtAssign(&o_y, m);
    q_y = NULL;
    m = emlrtCreateString1R2022a(&b_st, '1');
    emlrtAssign(&q_y, m);
    c_st.site = &ib_emlrtRSI;
    error(&c_st, b_y, getString(&c_st, b_message(&c_st, e_y, j_y, o_y, q_y,
            &f_emlrtMCI), &f_emlrtMCI), &f_emlrtMCI);
  }

  b_st.site = &m_emlrtRSI;
  p = true;
  if (b_value <= 0.0) {
    p = false;
  }

  if (!p) {
    for (i = 0; i < 23; i++) {
      g_u[i] = h_u[i];
    }

    d_y = NULL;
    m = emlrtCreateCharArray(2, &iv3[0]);
    emlrtInitCharArrayR2013a(&b_st, 23, m, &g_u[0]);
    emlrtAssign(&d_y, m);
    for (i = 0; i < 48; i++) {
      k_u[i] = m_u[i];
    }

    i_y = NULL;
    m = emlrtCreateCharArray(2, &iv7[0]);
    emlrtInitCharArrayR2013a(&b_st, 48, m, &k_u[0]);
    emlrtAssign(&i_y, m);
    for (i = 0; i < 23; i++) {
      g_u[i] = i_u[i];
    }

    n_y = NULL;
    m = emlrtCreateCharArray(2, &iv11[0]);
    emlrtInitCharArrayR2013a(&b_st, 23, m, &g_u[0]);
    emlrtAssign(&n_y, m);
    c_st.site = &gb_emlrtRSI;
    error(&c_st, d_y, getString(&c_st, message(&c_st, i_y, n_y, &emlrtMCI),
           &emlrtMCI), &emlrtMCI);
  }

  b_st.site = &m_emlrtRSI;
  p = true;
  if ((!!muDoubleScalarIsInf(b_value)) || (!!muDoubleScalarIsNaN(b_value))) {
    p = false;
  }

  if (!p) {
    for (i = 0; i < 21; i++) {
      j_u[i] = l_u[i];
    }

    g_y = NULL;
    m = emlrtCreateCharArray(2, &iv6[0]);
    emlrtInitCharArrayR2013a(&b_st, 21, m, &j_u[0]);
    emlrtAssign(&g_y, m);
    for (i = 0; i < 46; i++) {
      n_u[i] = q_u[i];
    }

    m_y = NULL;
    m = emlrtCreateCharArray(2, &iv10[0]);
    emlrtInitCharArrayR2013a(&b_st, 46, m, &n_u[0]);
    emlrtAssign(&m_y, m);
    for (i = 0; i < 23; i++) {
      g_u[i] = i_u[i];
    }

    r_y = NULL;
    m = emlrtCreateCharArray(2, &iv14[0]);
    emlrtInitCharArrayR2013a(&b_st, 23, m, &g_u[0]);
    emlrtAssign(&r_y, m);
    c_st.site = &fb_emlrtRSI;
    error(&c_st, g_y, getString(&c_st, message(&c_st, m_y, r_y, &b_emlrtMCI),
           &b_emlrtMCI), &b_emlrtMCI);
  }

  b_st.site = &m_emlrtRSI;
  p = true;
  if (muDoubleScalarIsNaN(b_value)) {
    p = false;
  }

  if (!p) {
    for (i = 0; i < 21; i++) {
      j_u[i] = o_u[i];
    }

    l_y = NULL;
    m = emlrtCreateCharArray(2, &iv9[0]);
    emlrtInitCharArrayR2013a(&b_st, 21, m, &j_u[0]);
    emlrtAssign(&l_y, m);
    for (i = 0; i < 46; i++) {
      n_u[i] = s_u[i];
    }

    p_y = NULL;
    m = emlrtCreateCharArray(2, &iv13[0]);
    emlrtInitCharArrayR2013a(&b_st, 46, m, &n_u[0]);
    emlrtAssign(&p_y, m);
    for (i = 0; i < 23; i++) {
      g_u[i] = i_u[i];
    }

    s_y = NULL;
    m = emlrtCreateCharArray(2, &iv15[0]);
    emlrtInitCharArrayR2013a(&b_st, 23, m, &g_u[0]);
    emlrtAssign(&s_y, m);
    c_st.site = &eb_emlrtRSI;
    error(&c_st, l_y, getString(&c_st, message(&c_st, p_y, s_y, &c_emlrtMCI),
           &c_emlrtMCI), &c_emlrtMCI);
  }

  obj->NormalizedLoopBandwidth = b_value;
}

static void CarrierSynchronizer_CalculateLoopGains(comm_CarrierSynchronizer *obj)
{
  static char_T b_b[4] = { 'A', 'u', 't', 'o' };

  real_T PhaseRecoveryGain;
  real_T PhaseRecoveryLoopBandwidth;
  real_T d;
  int32_T ret;
  char_T a[4];
  char_T b[4];
  PhaseRecoveryLoopBandwidth = obj->NormalizedLoopBandwidth *
    obj->SamplesPerSymbol;
  PhaseRecoveryGain = obj->SamplesPerSymbol;
  PhaseRecoveryLoopBandwidth /= (obj->DampingFactor + 0.25 / obj->DampingFactor)
    * obj->SamplesPerSymbol;
  d = (2.0 * obj->DampingFactor * PhaseRecoveryLoopBandwidth + 1.0) +
    PhaseRecoveryLoopBandwidth * PhaseRecoveryLoopBandwidth;
  obj->pProportionalGain = 4.0 * obj->DampingFactor * PhaseRecoveryLoopBandwidth
    / d / (2.0 * PhaseRecoveryGain);
  obj->pIntegratorGain = 4.0 * PhaseRecoveryLoopBandwidth *
    PhaseRecoveryLoopBandwidth / d / (2.0 * PhaseRecoveryGain);
  for (ret = 0; ret < 4; ret++) {
    a[ret] = obj->ModulationPhaseOffset[ret];
  }

  for (ret = 0; ret < 4; ret++) {
    b[ret] = b_b[ret];
  }

  ret = memcmp(&a[0], &b[0], 4);
  if (ret == 0) {
    obj->pActualPhaseOffset = 0.0;
  } else {
    obj->pActualPhaseOffset = obj->CustomPhaseOffset - 0.78539816339744828;
  }
}

static void SystemCore_checkTunablePropChange(const emlrtStack *sp,
  comm_CarrierSynchronizer *obj)
{
  static const int32_T iv[2] = { 1, 44 };

  static const int32_T iv1[2] = { 1, 44 };

  static char_T b_u[44] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 's', 'y', 's',
    't', 'e', 'm', ':', 'i', 'n', 'v', 'a', 'l', 'i', 'd', 'T', 'u', 'n', 'a',
    'b', 'l', 'e', 'M', 'o', 'd', 'A', 'c', 'c', 'e', 's', 's', 'C', 'o', 'd',
    'e', 'g', 'e', 'n' };

  emlrtStack st;
  const mxArray *b_y;
  const mxArray *m;
  const mxArray *y;
  int32_T i;
  char_T u[44];
  st.prev = sp;
  st.tls = sp->tls;
  if (obj->TunablePropsChanged) {
    for (i = 0; i < 44; i++) {
      u[i] = b_u[i];
    }

    y = NULL;
    m = emlrtCreateCharArray(2, &iv[0]);
    emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 44, m, &u[0]);
    emlrtAssign(&y, m);
    for (i = 0; i < 44; i++) {
      u[i] = b_u[i];
    }

    b_y = NULL;
    m = emlrtCreateCharArray(2, &iv1[0]);
    emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 44, m, &u[0]);
    emlrtAssign(&b_y, m);
    st.site = &e_emlrtRSI;
    error(&st, y, getString(&st, c_message(&st, b_y, &g_emlrtMCI), &g_emlrtMCI),
          &g_emlrtMCI);
  }
}

static void mw__internal__call__reset(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const emlrtStack *sp, real_T SamplesPerSymbol, real_T
  DampingFactor, real_T NormalizedLoopBandwidth)
{
  static const int32_T iv[2] = { 1, 45 };

  static const int32_T iv1[2] = { 1, 44 };

  static const int32_T iv2[2] = { 1, 45 };

  static const int32_T iv3[2] = { 1, 44 };

  static const int32_T iv4[2] = { 1, 5 };

  static char_T b_u[45] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 's', 'y', 's',
    't', 'e', 'm', ':', 'm', 'e', 't', 'h', 'o', 'd', 'C', 'a', 'l', 'l', 'e',
    'd', 'W', 'h', 'e', 'n', 'R', 'e', 'l', 'e', 'a', 's', 'e', 'd', 'C', 'o',
    'd', 'e', 'g', 'e', 'n' };

  static char_T d_u[44] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 's', 'y', 's',
    't', 'e', 'm', ':', 'i', 'n', 'v', 'a', 'l', 'i', 'd', 'T', 'u', 'n', 'a',
    'b', 'l', 'e', 'M', 'o', 'd', 'A', 'c', 'c', 'e', 's', 's', 'C', 'o', 'd',
    'e', 'g', 'e', 'n' };

  static char_T f_u[5] = { 'r', 'e', 's', 'e', 't' };

  static char_T correctedValue[4] = { 'A', 'u', 't', 'o' };

  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *m;
  const mxArray *y;
  int32_T i;
  char_T u[45];
  char_T c_u[44];
  char_T e_u[5];
  boolean_T tunablePropChangedBeforeResetImpl;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  if (!moduleInstance->sysobj_not_empty) {
    st.site = &f_emlrtRSI;
    CarrierSynchronizer_CarrierSynchronizer(&moduleInstance->sysobj);
    moduleInstance->sysobj_not_empty = true;
    st.site = &g_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    c_st.site = &d_emlrtRSI;
    tunablePropChangedBeforeResetImpl = (moduleInstance->sysobj.isInitialized ==
      1);
    if (tunablePropChangedBeforeResetImpl) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    for (i = 0; i < 4; i++) {
      moduleInstance->sysobj.ModulationPhaseOffset[i] = correctedValue[i];
    }

    st.site = &h_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    tunablePropChangedBeforeResetImpl = (moduleInstance->sysobj.isInitialized ==
      1);
    if (tunablePropChangedBeforeResetImpl) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    st.site = &h_emlrtRSI;
    CarrierSynchronizer_set_SamplesPerSymbol(&st, &moduleInstance->sysobj,
      SamplesPerSymbol);
    st.site = &i_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    tunablePropChangedBeforeResetImpl = (moduleInstance->sysobj.isInitialized ==
      1);
    if (tunablePropChangedBeforeResetImpl) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    st.site = &i_emlrtRSI;
    CarrierSynchronizer_set_DampingFactor(&st, &moduleInstance->sysobj,
      DampingFactor);
    st.site = &j_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    tunablePropChangedBeforeResetImpl = (moduleInstance->sysobj.isInitialized ==
      1);
    if (tunablePropChangedBeforeResetImpl) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    st.site = &j_emlrtRSI;
    CarrierSynchronizer_set_NormalizedLoopBandwidth(&st, &moduleInstance->sysobj,
      NormalizedLoopBandwidth);
  }

  st.site = &q_emlrtRSI;
  if (moduleInstance->sysobj.isInitialized == 2) {
    for (i = 0; i < 45; i++) {
      u[i] = b_u[i];
    }

    y = NULL;
    m = emlrtCreateCharArray(2, &iv[0]);
    emlrtInitCharArrayR2013a(&st, 45, m, &u[0]);
    emlrtAssign(&y, m);
    for (i = 0; i < 45; i++) {
      u[i] = b_u[i];
    }

    c_y = NULL;
    m = emlrtCreateCharArray(2, &iv2[0]);
    emlrtInitCharArrayR2013a(&st, 45, m, &u[0]);
    emlrtAssign(&c_y, m);
    for (i = 0; i < 5; i++) {
      e_u[i] = f_u[i];
    }

    e_y = NULL;
    m = emlrtCreateCharArray(2, &iv4[0]);
    emlrtInitCharArrayR2013a(&st, 5, m, &e_u[0]);
    emlrtAssign(&e_y, m);
    b_st.site = &e_emlrtRSI;
    error(&b_st, y, getString(&b_st, message(&b_st, c_y, e_y, &g_emlrtMCI),
           &g_emlrtMCI), &g_emlrtMCI);
  }

  tunablePropChangedBeforeResetImpl = moduleInstance->sysobj.TunablePropsChanged;
  if (moduleInstance->sysobj.isInitialized == 1) {
    b_st.site = &e_emlrtRSI;
    moduleInstance->sysobj.pLoopFilterState = 0.0;
    moduleInstance->sysobj.pIntegFilterState = 0.0;
    moduleInstance->sysobj.pDDSPreviousInput = 0.0;
    moduleInstance->sysobj.pPhase = 0.0;
    moduleInstance->sysobj.pPreviousSample.re = 0.0;
    moduleInstance->sysobj.pPreviousSample.im = 0.0;
  }

  if ((int32_T)tunablePropChangedBeforeResetImpl != (int32_T)
      moduleInstance->sysobj.TunablePropsChanged) {
    for (i = 0; i < 44; i++) {
      c_u[i] = d_u[i];
    }

    b_y = NULL;
    m = emlrtCreateCharArray(2, &iv1[0]);
    emlrtInitCharArrayR2013a(&st, 44, m, &c_u[0]);
    emlrtAssign(&b_y, m);
    for (i = 0; i < 44; i++) {
      c_u[i] = d_u[i];
    }

    d_y = NULL;
    m = emlrtCreateCharArray(2, &iv3[0]);
    emlrtInitCharArrayR2013a(&st, 44, m, &c_u[0]);
    emlrtAssign(&d_y, m);
    b_st.site = &e_emlrtRSI;
    error(&b_st, b_y, getString(&b_st, c_message(&b_st, d_y, &g_emlrtMCI),
           &g_emlrtMCI), &g_emlrtMCI);
  }
}

static void mw__internal__call__step(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const emlrtStack *sp, real_T SamplesPerSymbol, real_T
  DampingFactor, real_T NormalizedLoopBandwidth, coder_array_creal_T_2D *u0,
  coder_array_creal_T_2D *b_y0)
{
  coder_array_creal_T b_u0;
  coder_array_creal_T c_y0;
  emlrtStack st;
  int32_T i;
  int32_T u0_idx_0;
  st.prev = sp;
  st.tls = sp->tls;
  u0_idx_0 = u0->size[0];
  b_u0.vector.data = u0->vector.data;
  b_u0.vector.numel = u0->vector.numel;
  b_u0.vector.allocated = u0->vector.allocated;
  b_u0.vector.owner = false;
  b_u0.size[0] = u0_idx_0;
  array_creal_T_Constructor(&c_y0);
  st.site = &db_emlrtRSI;
  mw__internal__system___fcn(moduleInstance, &st, SamplesPerSymbol,
    DampingFactor, NormalizedLoopBandwidth, &b_u0, &c_y0);
  emlrtDimSizeGeqCheckR2012b(6175, c_y0.size[0], &emlrtECI, (void *)sp);
  array_creal_T_2D_SetSize(b_y0, c_y0.size[0], 1);
  u0_idx_0 = c_y0.size[0];
  for (i = 0; i < u0_idx_0; i++) {
    b_y0->vector.data[i].re = c_y0.vector.data[i].re;
    b_y0->vector.data[i].im = c_y0.vector.data[i].im;
  }

  array_creal_T_Destructor(&c_y0);
}

static void mw__internal__system___fcn(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const emlrtStack *sp, real_T varargin_3, real_T varargin_4,
  real_T varargin_5, coder_array_creal_T *varargin_7, coder_array_creal_T
  *varargout_1)
{
  static const int32_T iv[2] = { 1, 51 };

  static const int32_T iv1[2] = { 1, 51 };

  static const int32_T iv2[2] = { 1, 5 };

  static char_T b_u[51] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 's', 'y', 's',
    't', 'e', 'm', ':', 'm', 'e', 't', 'h', 'o', 'd', 'C', 'a', 'l', 'l', 'e',
    'd', 'W', 'h', 'e', 'n', 'L', 'o', 'c', 'k', 'e', 'd', 'R', 'e', 'l', 'e',
    'a', 's', 'e', 'd', 'C', 'o', 'd', 'e', 'g', 'e', 'n' };

  static char_T d_u[5] = { 's', 'e', 't', 'u', 'p' };

  static char_T correctedValue[4] = { 'A', 'u', 't', 'o' };

  cell_wrap varSizes[1];
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack st;
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *m;
  const mxArray *y;
  real_T PhaseRecoveryGain;
  real_T PhaseRecoveryLoopBandwidth;
  real_T d;
  int32_T ret;
  int16_T inSize[8];
  char_T u[51];
  char_T c_u[5];
  char_T a[4];
  char_T b[4];
  boolean_T exitg1;
  boolean_T flag;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  if (!moduleInstance->sysobj_not_empty) {
    st.site = &f_emlrtRSI;
    CarrierSynchronizer_CarrierSynchronizer(&moduleInstance->sysobj);
    moduleInstance->sysobj_not_empty = true;
    st.site = &g_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    c_st.site = &d_emlrtRSI;
    flag = (moduleInstance->sysobj.isInitialized == 1);
    if (flag) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    for (ret = 0; ret < 4; ret++) {
      moduleInstance->sysobj.ModulationPhaseOffset[ret] = correctedValue[ret];
    }

    st.site = &h_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    flag = (moduleInstance->sysobj.isInitialized == 1);
    if (flag) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    st.site = &h_emlrtRSI;
    CarrierSynchronizer_set_SamplesPerSymbol(&st, &moduleInstance->sysobj,
      varargin_3);
    st.site = &i_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    flag = (moduleInstance->sysobj.isInitialized == 1);
    if (flag) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    st.site = &i_emlrtRSI;
    CarrierSynchronizer_set_DampingFactor(&st, &moduleInstance->sysobj,
      varargin_4);
    st.site = &j_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    flag = (moduleInstance->sysobj.isInitialized == 1);
    if (flag) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    st.site = &j_emlrtRSI;
    CarrierSynchronizer_set_NormalizedLoopBandwidth(&st, &moduleInstance->sysobj,
      varargin_5);
  }

  if (moduleInstance->sysobj.SamplesPerSymbol != varargin_3) {
    st.site = &r_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    flag = (moduleInstance->sysobj.isInitialized == 1);
    if (flag) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    st.site = &r_emlrtRSI;
    CarrierSynchronizer_set_SamplesPerSymbol(&st, &moduleInstance->sysobj,
      varargin_3);
  }

  if (moduleInstance->sysobj.DampingFactor != varargin_4) {
    st.site = &s_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    flag = (moduleInstance->sysobj.isInitialized == 1);
    if (flag) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    st.site = &s_emlrtRSI;
    CarrierSynchronizer_set_DampingFactor(&st, &moduleInstance->sysobj,
      varargin_4);
  }

  if (moduleInstance->sysobj.NormalizedLoopBandwidth != varargin_5) {
    st.site = &t_emlrtRSI;
    b_st.site = &d_emlrtRSI;
    flag = (moduleInstance->sysobj.isInitialized == 1);
    if (flag) {
      moduleInstance->sysobj.TunablePropsChanged = true;
    }

    st.site = &t_emlrtRSI;
    CarrierSynchronizer_set_NormalizedLoopBandwidth(&st, &moduleInstance->sysobj,
      varargin_5);
  }

  st.site = &u_emlrtRSI;
  if (moduleInstance->sysobj.isInitialized != 1) {
    b_st.site = &e_emlrtRSI;
    c_st.site = &e_emlrtRSI;
    if (moduleInstance->sysobj.isInitialized != 0) {
      for (ret = 0; ret < 51; ret++) {
        u[ret] = b_u[ret];
      }

      y = NULL;
      m = emlrtCreateCharArray(2, &iv[0]);
      emlrtInitCharArrayR2013a(&c_st, 51, m, &u[0]);
      emlrtAssign(&y, m);
      for (ret = 0; ret < 51; ret++) {
        u[ret] = b_u[ret];
      }

      b_y = NULL;
      m = emlrtCreateCharArray(2, &iv1[0]);
      emlrtInitCharArrayR2013a(&c_st, 51, m, &u[0]);
      emlrtAssign(&b_y, m);
      for (ret = 0; ret < 5; ret++) {
        c_u[ret] = d_u[ret];
      }

      c_y = NULL;
      m = emlrtCreateCharArray(2, &iv2[0]);
      emlrtInitCharArrayR2013a(&c_st, 5, m, &c_u[0]);
      emlrtAssign(&c_y, m);
      d_st.site = &e_emlrtRSI;
      error(&d_st, y, getString(&d_st, message(&d_st, b_y, c_y, &g_emlrtMCI),
             &g_emlrtMCI), &g_emlrtMCI);
    }

    moduleInstance->sysobj.isInitialized = 1;
    d_st.site = &e_emlrtRSI;
    varSizes[0].f1[0] = (uint32_T)varargin_7->size[0];
    varSizes[0].f1[1] = 1U;
    for (ret = 0; ret < 6; ret++) {
      varSizes[0].f1[ret + 2] = 1U;
    }

    moduleInstance->sysobj.inputVarSize[0] = varSizes[0];
    d_st.site = &e_emlrtRSI;
    moduleInstance->sysobj.pPhase = 0.0;
    moduleInstance->sysobj.pPreviousSample.re = 0.0;
    moduleInstance->sysobj.pPreviousSample.im = 0.0;
    e_st.site = &p_emlrtRSI;
    CarrierSynchronizer_CalculateLoopGains(&moduleInstance->sysobj);
    moduleInstance->sysobj.pDigitalSynthesizerGain = -1.0;
    d_st.site = &e_emlrtRSI;
    SystemCore_checkTunablePropChange(&d_st, &moduleInstance->sysobj);
    moduleInstance->sysobj.TunablePropsChanged = false;
    c_st.site = &e_emlrtRSI;
    moduleInstance->sysobj.pLoopFilterState = 0.0;
    moduleInstance->sysobj.pIntegFilterState = 0.0;
    moduleInstance->sysobj.pDDSPreviousInput = 0.0;
    moduleInstance->sysobj.pPhase = 0.0;
    moduleInstance->sysobj.pPreviousSample.re = 0.0;
    moduleInstance->sysobj.pPreviousSample.im = 0.0;
  }

  b_st.site = &e_emlrtRSI;
  if (!moduleInstance->sysobj.CacheInputSizes) {
    moduleInstance->sysobj.CacheInputSizes = true;
    c_st.site = &e_emlrtRSI;
    varSizes[0].f1[0] = (uint32_T)varargin_7->size[0];
    varSizes[0].f1[1] = 1U;
    for (ret = 0; ret < 6; ret++) {
      varSizes[0].f1[ret + 2] = 1U;
    }

    moduleInstance->sysobj.inputVarSize[0] = varSizes[0];
  }

  b_st.site = &e_emlrtRSI;
  if (moduleInstance->sysobj.TunablePropsChanged) {
    moduleInstance->sysobj.TunablePropsChanged = false;
    c_st.site = &e_emlrtRSI;
    d_st.site = &v_emlrtRSI;
    PhaseRecoveryLoopBandwidth = moduleInstance->sysobj.NormalizedLoopBandwidth *
      moduleInstance->sysobj.SamplesPerSymbol;
    PhaseRecoveryGain = moduleInstance->sysobj.SamplesPerSymbol;
    PhaseRecoveryLoopBandwidth /= (moduleInstance->sysobj.DampingFactor + 0.25 /
      moduleInstance->sysobj.DampingFactor) *
      moduleInstance->sysobj.SamplesPerSymbol;
    d = (2.0 * moduleInstance->sysobj.DampingFactor * PhaseRecoveryLoopBandwidth
         + 1.0) + PhaseRecoveryLoopBandwidth * PhaseRecoveryLoopBandwidth;
    moduleInstance->sysobj.pProportionalGain = 4.0 *
      moduleInstance->sysobj.DampingFactor * PhaseRecoveryLoopBandwidth / d /
      (2.0 * PhaseRecoveryGain);
    moduleInstance->sysobj.pIntegratorGain = 4.0 * PhaseRecoveryLoopBandwidth *
      PhaseRecoveryLoopBandwidth / d / (2.0 * PhaseRecoveryGain);
    for (ret = 0; ret < 4; ret++) {
      a[ret] = moduleInstance->sysobj.ModulationPhaseOffset[ret];
    }

    for (ret = 0; ret < 4; ret++) {
      b[ret] = correctedValue[ret];
    }

    ret = memcmp(&a[0], &b[0], 4);
    if (ret == 0) {
      moduleInstance->sysobj.pActualPhaseOffset = 0.0;
    } else {
      moduleInstance->sysobj.pActualPhaseOffset =
        moduleInstance->sysobj.CustomPhaseOffset - 0.78539816339744828;
    }
  }

  b_st.site = &e_emlrtRSI;
  inSize[0] = (int16_T)varargin_7->size[0];
  inSize[1] = 1;
  for (ret = 0; ret < 6; ret++) {
    inSize[ret + 2] = 1;
  }

  ret = 0;
  exitg1 = false;
  while ((!exitg1) && (ret < 8)) {
    if (moduleInstance->sysobj.inputVarSize[0].f1[ret] != (uint32_T)inSize[ret])
    {
      for (ret = 0; ret < 8; ret++) {
        moduleInstance->sysobj.inputVarSize[0].f1[ret] = (uint32_T)inSize[ret];
      }

      exitg1 = true;
    } else {
      ret++;
    }
  }

  b_st.site = &e_emlrtRSI;
  CarrierSynchronizer_stepImpl(&b_st, &moduleInstance->sysobj, varargin_7,
    varargout_1);
  b_st.site = &e_emlrtRSI;
  SystemCore_checkTunablePropChange(&b_st, &moduleInstance->sysobj);
}

static void CarrierSynchronizer_stepImpl(const emlrtStack *sp,
  comm_CarrierSynchronizer *obj, coder_array_creal_T *input, coder_array_creal_T
  *output)
{
  emlrtStack st;
  creal_T previousSample;
  real_T phaseEstimate_data[6175];
  real_T DDSPreviousInp;
  real_T digitalSynthesizerGain;
  real_T integFiltState;
  real_T integratorGain;
  real_T loopFiltState;
  real_T phase;
  real_T phaseOffset;
  real_T proportionalGain;
  int32_T phaseEstimate_size[1];
  st.prev = sp;
  st.tls = sp->tls;
  loopFiltState = obj->pLoopFilterState;
  integFiltState = obj->pIntegFilterState;
  DDSPreviousInp = obj->pDDSPreviousInput;
  previousSample.re = obj->pPreviousSample.re;
  previousSample.im = obj->pPreviousSample.im;
  phase = obj->pPhase;
  digitalSynthesizerGain = obj->pDigitalSynthesizerGain;
  integratorGain = obj->pIntegratorGain;
  proportionalGain = obj->pProportionalGain;
  phaseOffset = obj->pActualPhaseOffset;
  st.site = &w_emlrtRSI;
  carrierSyncCore(&st, input, &loopFiltState, &integFiltState, &DDSPreviousInp,
                  &previousSample, &phase, &digitalSynthesizerGain,
                  &integratorGain, &proportionalGain, &phaseOffset, output,
                  phaseEstimate_data, phaseEstimate_size);
  obj->pLoopFilterState = loopFiltState;
  obj->pIntegFilterState = integFiltState;
  obj->pDDSPreviousInput = DDSPreviousInp;
  obj->pPreviousSample.re = previousSample.re;
  obj->pPreviousSample.im = previousSample.im;
  obj->pPhase = phase;
  obj->pDigitalSynthesizerGain = digitalSynthesizerGain;
  obj->pIntegratorGain = integratorGain;
  obj->pProportionalGain = proportionalGain;
  obj->pActualPhaseOffset = phaseOffset;
}

static creal_T expFunc(real_T in)
{
  creal_T val;
  real_T cosVal;
  real_T sinVal;
  sinVal = sin(in);
  cosVal = cos(in);
  val.re = cosVal;
  val.im = sinVal;
  return val;
}

static const mxArray *emlrt_marshallOut(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const emlrtStack *sp)
{
  static const int32_T iv2[2] = { 1, 8 };

  static const int32_T iv3[2] = { 1, 4 };

  static const int32_T iv[1] = { 2 };

  static const char_T *sv[18] = { "isInitialized", "TunablePropsChanged",
    "inputVarSize", "CacheInputSizes", "ModulationPhaseOffset",
    "CustomPhaseOffset", "SamplesPerSymbol", "DampingFactor",
    "NormalizedLoopBandwidth", "pProportionalGain", "pIntegratorGain",
    "pDigitalSynthesizerGain", "pPhase", "pPreviousSample", "pActualPhaseOffset",
    "pLoopFilterState", "pIntegFilterState", "pDDSPreviousInput" };

  static const char_T *sv1[1] = { "f1" };

  cell_wrap c_u[1];
  const mxArray *b_y;
  const mxArray *c_y;
  const mxArray *d_y;
  const mxArray *e_y;
  const mxArray *f_y;
  const mxArray *g_y;
  const mxArray *h_y;
  const mxArray *i_y;
  const mxArray *j_y;
  const mxArray *k_y;
  const mxArray *l_y;
  const mxArray *m;
  const mxArray *m_y;
  const mxArray *n_y;
  const mxArray *o_y;
  const mxArray *p_y;
  const mxArray *q_y;
  const mxArray *r_y;
  const mxArray *s_y;
  const mxArray *t_y;
  const mxArray *u_y;
  const mxArray *v_y;
  const mxArray *w_y;
  const mxArray *y;
  creal_T *r;
  real_T f_u;
  real_T u_im;
  int32_T iv1[1];
  int32_T i;
  int32_T u;
  int32_T *pData;
  uint32_T d_u[8];
  uint32_T *b_pData;
  char_T e_u[4];
  boolean_T b_u;
  y = NULL;
  emlrtAssign(&y, emlrtCreateCellMatrix(3, 1));
  b_y = NULL;
  m = emlrtCreateNumericArray(1, (const void *)&iv[0], mxINT32_CLASS, mxREAL);
  pData = (int32_T *)emlrtMxGetData(m);
  u = 0;
  for (i = 0; i < 2; i++) {
    pData[u] = (*moduleInstance->v)[i];
    u++;
  }

  emlrtAssign(&b_y, m);
  emlrtSetCell(y, 0, b_y);
  c_y = NULL;
  emlrtAssign(&c_y, emlrtCreateStructMatrix(1, 1, 18, (const char_T **)&sv[0]));
  u = moduleInstance->sysobj.isInitialized;
  d_y = NULL;
  m = emlrtCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
  *(int32_T *)emlrtMxGetData(m) = u;
  emlrtAssign(&d_y, m);
  emlrtSetFieldR2017b(c_y, 0, "isInitialized", d_y, 0);
  b_u = moduleInstance->sysobj.TunablePropsChanged;
  e_y = NULL;
  m = emlrtCreateLogicalScalar(b_u);
  emlrtAssign(&e_y, m);
  emlrtSetFieldR2017b(c_y, 0, "TunablePropsChanged", e_y, 1);
  c_u[0] = moduleInstance->sysobj.inputVarSize[0];
  f_y = NULL;
  iv1[0] = 1;
  emlrtAssign(&f_y, emlrtCreateStructArray(1, &iv1[0], 1, (const char_T **)&sv1
    [0]));
  for (u = 0; u < 8; u++) {
    d_u[u] = c_u[0].f1[u];
  }

  g_y = NULL;
  m = emlrtCreateNumericArray(2, (const void *)&iv2[0], mxUINT32_CLASS, mxREAL);
  b_pData = (uint32_T *)emlrtMxGetData(m);
  u = 0;
  for (i = 0; i < 8; i++) {
    b_pData[u] = d_u[i];
    u++;
  }

  emlrtAssign(&g_y, m);
  emlrtSetFieldR2017b(f_y, 0, "f1", g_y, 0);
  emlrtSetFieldR2017b(c_y, 0, "inputVarSize", f_y, 2);
  b_u = moduleInstance->sysobj.CacheInputSizes;
  h_y = NULL;
  m = emlrtCreateLogicalScalar(b_u);
  emlrtAssign(&h_y, m);
  emlrtSetFieldR2017b(c_y, 0, "CacheInputSizes", h_y, 3);
  for (u = 0; u < 4; u++) {
    e_u[u] = moduleInstance->sysobj.ModulationPhaseOffset[u];
  }

  i_y = NULL;
  m = emlrtCreateCharArray(2, &iv3[0]);
  emlrtInitCharArrayR2013a((emlrtConstCTX)sp, 4, m, &e_u[0]);
  emlrtAssign(&i_y, m);
  emlrtSetFieldR2017b(c_y, 0, "ModulationPhaseOffset", i_y, 4);
  f_u = moduleInstance->sysobj.CustomPhaseOffset;
  j_y = NULL;
  m = emlrtCreateDoubleScalar(f_u);
  emlrtAssign(&j_y, m);
  emlrtSetFieldR2017b(c_y, 0, "CustomPhaseOffset", j_y, 5);
  f_u = moduleInstance->sysobj.SamplesPerSymbol;
  k_y = NULL;
  m = emlrtCreateDoubleScalar(f_u);
  emlrtAssign(&k_y, m);
  emlrtSetFieldR2017b(c_y, 0, "SamplesPerSymbol", k_y, 6);
  f_u = moduleInstance->sysobj.DampingFactor;
  l_y = NULL;
  m = emlrtCreateDoubleScalar(f_u);
  emlrtAssign(&l_y, m);
  emlrtSetFieldR2017b(c_y, 0, "DampingFactor", l_y, 7);
  f_u = moduleInstance->sysobj.NormalizedLoopBandwidth;
  m_y = NULL;
  m = emlrtCreateDoubleScalar(f_u);
  emlrtAssign(&m_y, m);
  emlrtSetFieldR2017b(c_y, 0, "NormalizedLoopBandwidth", m_y, 8);
  f_u = moduleInstance->sysobj.pProportionalGain;
  n_y = NULL;
  m = emlrtCreateDoubleScalar(f_u);
  emlrtAssign(&n_y, m);
  emlrtSetFieldR2017b(c_y, 0, "pProportionalGain", n_y, 9);
  f_u = moduleInstance->sysobj.pIntegratorGain;
  o_y = NULL;
  m = emlrtCreateDoubleScalar(f_u);
  emlrtAssign(&o_y, m);
  emlrtSetFieldR2017b(c_y, 0, "pIntegratorGain", o_y, 10);
  f_u = moduleInstance->sysobj.pDigitalSynthesizerGain;
  p_y = NULL;
  m = emlrtCreateDoubleScalar(f_u);
  emlrtAssign(&p_y, m);
  emlrtSetFieldR2017b(c_y, 0, "pDigitalSynthesizerGain", p_y, 11);
  f_u = moduleInstance->sysobj.pPhase;
  q_y = NULL;
  m = emlrtCreateDoubleScalar(f_u);
  emlrtAssign(&q_y, m);
  emlrtSetFieldR2017b(c_y, 0, "pPhase", q_y, 12);
  f_u = moduleInstance->sysobj.pPreviousSample.re;
  u_im = moduleInstance->sysobj.pPreviousSample.im;
  r_y = NULL;
  m = emlrtCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxCOMPLEX);
  r = (creal_T *)emlrtMxGetData(m);
  r->re = f_u;
  r->im = u_im;
  emlrtFreeImagIfZero(m);
  emlrtAssign(&r_y, m);
  emlrtSetFieldR2017b(c_y, 0, "pPreviousSample", r_y, 13);
  f_u = moduleInstance->sysobj.pActualPhaseOffset;
  s_y = NULL;
  m = emlrtCreateDoubleScalar(f_u);
  emlrtAssign(&s_y, m);
  emlrtSetFieldR2017b(c_y, 0, "pActualPhaseOffset", s_y, 14);
  f_u = moduleInstance->sysobj.pLoopFilterState;
  t_y = NULL;
  m = emlrtCreateDoubleScalar(f_u);
  emlrtAssign(&t_y, m);
  emlrtSetFieldR2017b(c_y, 0, "pLoopFilterState", t_y, 15);
  f_u = moduleInstance->sysobj.pIntegFilterState;
  u_y = NULL;
  m = emlrtCreateDoubleScalar(f_u);
  emlrtAssign(&u_y, m);
  emlrtSetFieldR2017b(c_y, 0, "pIntegFilterState", u_y, 16);
  f_u = moduleInstance->sysobj.pDDSPreviousInput;
  v_y = NULL;
  m = emlrtCreateDoubleScalar(f_u);
  emlrtAssign(&v_y, m);
  emlrtSetFieldR2017b(c_y, 0, "pDDSPreviousInput", v_y, 17);
  emlrtSetCell(y, 1, c_y);
  w_y = NULL;
  m = emlrtCreateLogicalScalar(moduleInstance->sysobj_not_empty);
  emlrtAssign(&w_y, m);
  emlrtSetCell(y, 2, w_y);
  return y;
}

static const mxArray *cgxe_mdl_get_sim_state
  (InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance)
{
  emlrtStack b_st = { NULL,            /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  const mxArray *st;
  b_st.tls = moduleInstance->emlrtRootTLSGlobal;
  st = NULL;
  emlrtAssign(&st, emlrt_marshallOut(moduleInstance, &b_st));
  return st;
}

static void emlrt_marshallIn(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const emlrtStack *sp, const mxArray *u)
{
  emlrtMsgIdentifier thisId;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  thisId.fIdentifier = "v";
  b_emlrt_marshallIn(sp, emlrtAlias(emlrtGetCell((void *)sp, &thisId, u, 0)),
                     "v", *moduleInstance->v);
  thisId.fIdentifier = "sysobj";
  d_emlrt_marshallIn(sp, emlrtAlias(emlrtGetCell((void *)sp, &thisId, u, 1)),
                     "sysobj", &moduleInstance->sysobj);
  thisId.fIdentifier = "sysobj_not_empty";
  moduleInstance->sysobj_not_empty = m_emlrt_marshallIn(sp, emlrtAlias
    (emlrtGetCell((void *)sp, &thisId, u, 2)), "sysobj_not_empty");
  emlrtDestroyArray(&u);
}

static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
  const char_T *identifier, int32_T y[2])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  c_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId, y);
  emlrtDestroyArray(&nullptr);
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, int32_T y[2])
{
  n_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
  const char_T *identifier, comm_CarrierSynchronizer *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  e_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId, y);
  emlrtDestroyArray(&nullptr);
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, comm_CarrierSynchronizer *y)
{
  static const int32_T dims = 0;
  static const char_T *fieldNames[18] = { "isInitialized", "TunablePropsChanged",
    "inputVarSize", "CacheInputSizes", "ModulationPhaseOffset",
    "CustomPhaseOffset", "SamplesPerSymbol", "DampingFactor",
    "NormalizedLoopBandwidth", "pProportionalGain", "pIntegratorGain",
    "pDigitalSynthesizerGain", "pPhase", "pPreviousSample", "pActualPhaseOffset",
    "pLoopFilterState", "pIntegFilterState", "pDDSPreviousInput" };

  emlrtMsgIdentifier thisId;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 18, (const char_T **)
    &fieldNames[0], 0U, (const void *)&dims);
  thisId.fIdentifier = "isInitialized";
  y->isInitialized = f_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 0, "isInitialized")), &thisId);
  thisId.fIdentifier = "TunablePropsChanged";
  y->TunablePropsChanged = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 1, "TunablePropsChanged")), &thisId);
  thisId.fIdentifier = "inputVarSize";
  h_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0,
    2, "inputVarSize")), &thisId, y->inputVarSize);
  thisId.fIdentifier = "CacheInputSizes";
  y->CacheInputSizes = g_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 3, "CacheInputSizes")), &thisId);
  thisId.fIdentifier = "ModulationPhaseOffset";
  j_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0,
    4, "ModulationPhaseOffset")), &thisId, y->ModulationPhaseOffset);
  thisId.fIdentifier = "CustomPhaseOffset";
  y->CustomPhaseOffset = k_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 5, "CustomPhaseOffset")), &thisId);
  thisId.fIdentifier = "SamplesPerSymbol";
  y->SamplesPerSymbol = k_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 6, "SamplesPerSymbol")), &thisId);
  thisId.fIdentifier = "DampingFactor";
  y->DampingFactor = k_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 7, "DampingFactor")), &thisId);
  thisId.fIdentifier = "NormalizedLoopBandwidth";
  y->NormalizedLoopBandwidth = k_emlrt_marshallIn(sp, emlrtAlias
    (emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 8, "NormalizedLoopBandwidth")),
    &thisId);
  thisId.fIdentifier = "pProportionalGain";
  y->pProportionalGain = k_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 9, "pProportionalGain")), &thisId);
  thisId.fIdentifier = "pIntegratorGain";
  y->pIntegratorGain = k_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 10, "pIntegratorGain")), &thisId);
  thisId.fIdentifier = "pDigitalSynthesizerGain";
  y->pDigitalSynthesizerGain = k_emlrt_marshallIn(sp, emlrtAlias
    (emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0, 11, "pDigitalSynthesizerGain")),
    &thisId);
  thisId.fIdentifier = "pPhase";
  y->pPhase = k_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 12, "pPhase")), &thisId);
  thisId.fIdentifier = "pPreviousSample";
  y->pPreviousSample = l_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 13, "pPreviousSample")), &thisId);
  thisId.fIdentifier = "pActualPhaseOffset";
  y->pActualPhaseOffset = k_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 14, "pActualPhaseOffset")), &thisId);
  thisId.fIdentifier = "pLoopFilterState";
  y->pLoopFilterState = k_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 15, "pLoopFilterState")), &thisId);
  thisId.fIdentifier = "pIntegFilterState";
  y->pIntegFilterState = k_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 16, "pIntegFilterState")), &thisId);
  thisId.fIdentifier = "pDDSPreviousInput";
  y->pDDSPreviousInput = k_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b
    ((emlrtConstCTX)sp, u, 0, 17, "pDDSPreviousInput")), &thisId);
  emlrtDestroyArray(&u);
}

static int32_T f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  int32_T y;
  y = o_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static boolean_T g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u,
  const emlrtMsgIdentifier *parentId)
{
  boolean_T y;
  y = p_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, cell_wrap y[1])
{
  static const int32_T dims[1] = { 1 };

  static const char_T *fieldNames[1] = { "f1" };

  emlrtMsgIdentifier thisId;
  thisId.fParent = parentId;
  thisId.bParentIsCell = false;
  emlrtCheckStructR2012b((emlrtConstCTX)sp, parentId, u, 1, (const char_T **)
    &fieldNames[0], 1U, (const void *)&dims[0]);
  thisId.fIdentifier = "f1";
  i_emlrt_marshallIn(sp, emlrtAlias(emlrtGetFieldR2017b((emlrtConstCTX)sp, u, 0,
    0, "f1")), &thisId, y[0].f1);
  emlrtDestroyArray(&u);
}

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, uint32_T y[8])
{
  q_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[4])
{
  r_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  real_T y;
  y = s_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static creal_T l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId)
{
  creal_T y;
  y = t_emlrt_marshallIn(sp, emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}

static boolean_T m_emlrt_marshallIn(const emlrtStack *sp, const mxArray *nullptr,
  const char_T *identifier)
{
  emlrtMsgIdentifier thisId;
  boolean_T y;
  thisId.fIdentifier = (const char_T *)identifier;
  thisId.fParent = NULL;
  thisId.bParentIsCell = false;
  y = g_emlrt_marshallIn(sp, emlrtAlias(nullptr), &thisId);
  emlrtDestroyArray(&nullptr);
  return y;
}

static void cgxe_mdl_set_sim_state(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance, const mxArray *st)
{
  emlrtStack b_st = { NULL,            /* site */
    NULL,                              /* tls */
    NULL                               /* prev */
  };

  b_st.tls = moduleInstance->emlrtRootTLSGlobal;
  emlrt_marshallIn(moduleInstance, &b_st, emlrtAlias(st));
  emlrtDestroyArray(&st);
}

static const mxArray *message(const emlrtStack *sp, const mxArray *m1, const
  mxArray *m2, emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  const mxArray *m;
  pArrays[0] = m1;
  pArrays[1] = m2;
  return emlrtCallMATLABR2012b((emlrtConstCTX)sp, 1, &m, 2, &pArrays[0],
    "message", true, location);
}

static const mxArray *getString(const emlrtStack *sp, const mxArray *m1,
  emlrtMCInfo *location)
{
  const mxArray *m;
  const mxArray *pArray;
  pArray = m1;
  return emlrtCallMATLABR2012b((emlrtConstCTX)sp, 1, &m, 1, &pArray, "getString",
    true, location);
}

static void error(const emlrtStack *sp, const mxArray *m, const mxArray *m1,
                  emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  pArrays[0] = m;
  pArrays[1] = m1;
  emlrtCallMATLABR2012b((emlrtConstCTX)sp, 0, NULL, 2, &pArrays[0], "error",
                        true, location);
}

static const mxArray *b_message(const emlrtStack *sp, const mxArray *m1, const
  mxArray *m2, const mxArray *m3, const mxArray *m4, emlrtMCInfo *location)
{
  const mxArray *pArrays[4];
  const mxArray *m;
  pArrays[0] = m1;
  pArrays[1] = m2;
  pArrays[2] = m3;
  pArrays[3] = m4;
  return emlrtCallMATLABR2012b((emlrtConstCTX)sp, 1, &m, 4, &pArrays[0],
    "message", true, location);
}

static const mxArray *c_message(const emlrtStack *sp, const mxArray *m1,
  emlrtMCInfo *location)
{
  const mxArray *m;
  const mxArray *pArray;
  pArray = m1;
  return emlrtCallMATLABR2012b((emlrtConstCTX)sp, 1, &m, 1, &pArray, "message",
    true, location);
}

static void n_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, int32_T ret[2])
{
  static const int32_T dims[1] = { 2 };

  int32_T (*r)[2];
  int32_T i;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "int32", false, 1U, (
    const void *)&dims[0]);
  r = (int32_T (*)[2])emlrtMxGetData(src);
  for (i = 0; i < 2; i++) {
    ret[i] = (*r)[i];
  }

  emlrtDestroyArray(&src);
}

static int32_T o_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  int32_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "int32", false, 0U, (
    const void *)&dims);
  ret = *(int32_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static boolean_T p_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  boolean_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "logical", false, 0U, (
    const void *)&dims);
  ret = *emlrtMxGetLogicals(src);
  emlrtDestroyArray(&src);
  return ret;
}

static void q_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, uint32_T ret[8])
{
  static const int32_T dims[2] = { 1, 8 };

  int32_T i;
  uint32_T (*r)[8];
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "uint32", false, 2U, (
    const void *)&dims[0]);
  r = (uint32_T (*)[8])emlrtMxGetData(src);
  for (i = 0; i < 8; i++) {
    ret[i] = (*r)[i];
  }

  emlrtDestroyArray(&src);
}

static void r_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[4])
{
  static const int32_T dims[2] = { 1, 4 };

  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "char", false, 2U, (
    const void *)&dims[0]);
  emlrtImportCharArrayR2015b((emlrtConstCTX)sp, src, &ret[0], 4);
  emlrtDestroyArray(&src);
}

static real_T s_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  real_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", false, 0U, (
    const void *)&dims);
  ret = *(real_T *)emlrtMxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}

static creal_T t_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src,
  const emlrtMsgIdentifier *msgId)
{
  static const int32_T dims = 0;
  creal_T ret;
  emlrtCheckBuiltInR2012b((emlrtConstCTX)sp, msgId, src, "double", true, 0U, (
    const void *)&dims);
  emlrtImportArrayR2015b((emlrtConstCTX)sp, src, &ret, 8, true);
  emlrtDestroyArray(&src);
  return ret;
}

static void carrierSyncCore(const emlrtStack *sp, coder_array_creal_T *input,
  real_T *loopFiltState, real_T *integFiltState, real_T *DDSPreviousInp, creal_T
  *previousSample, real_T *phase, real_T *digitalSynthesizerGain, real_T
  *integratorGain, real_T *proportionalGain, real_T *phaseOffset,
  coder_array_creal_T *output, real_T phaseEstimate_data[], int32_T
  phaseEstimate_size[1])
{
  coder_array_creal_T phaseCorrection;
  emlrtStack b_st;
  emlrtStack st;
  creal_T b;
  real_T DDSOut;
  real_T loopFiltOut;
  real_T phErr;
  int32_T i;
  int32_T loop_ub;
  array_creal_T_SetSize(output, input->size[0]);
  loop_ub = input->size[0];
  for (i = 0; i < loop_ub; i++) {
    output->vector.data[i].re = 0.0;
    output->vector.data[i].im = 0.0;
  }

  array_creal_T_Constructor(&phaseCorrection);
  array_creal_T_SetSize(&phaseCorrection, input->size[0]);
  loop_ub = input->size[0];
  for (i = 0; i < loop_ub; i++) {
    phaseCorrection.vector.data[i].re = 0.0;
    phaseCorrection.vector.data[i].im = 0.0;
  }

  i = input->size[0] - 1;
  for (loop_ub = 0; loop_ub <= i; loop_ub++) {
    st.site = &x_emlrtRSI;
    b_st.site = &cb_emlrtRSI;
    loopFiltOut = signDouble(previousSample->re);
    st.site = &y_emlrtRSI;
    b_st.site = &cb_emlrtRSI;
    DDSOut = signDouble(previousSample->im);
    phErr = loopFiltOut * previousSample->im - DDSOut * previousSample->re;
    st.site = &ab_emlrtRSI;
    b = expFunc(*phase);
    loopFiltOut = input->vector.data[emlrtDynamicBoundsCheckR2012b(loop_ub + 1,
      1, input->size[0], &d_emlrtBCI, (emlrtConstCTX)sp) - 1].re * b.re -
      input->vector.data[emlrtDynamicBoundsCheckR2012b(loop_ub + 1, 1,
      input->size[0], &d_emlrtBCI, (emlrtConstCTX)sp) - 1].im * b.im;
    DDSOut = input->vector.data[emlrtDynamicBoundsCheckR2012b(loop_ub + 1, 1,
      input->size[0], &d_emlrtBCI, (emlrtConstCTX)sp) - 1].re * b.im +
      input->vector.data[emlrtDynamicBoundsCheckR2012b(loop_ub + 1, 1,
      input->size[0], &d_emlrtBCI, (emlrtConstCTX)sp) - 1].im * b.re;
    output->vector.data[emlrtDynamicBoundsCheckR2012b(loop_ub + 1, 1,
      output->size[0], &emlrtBCI, (emlrtConstCTX)sp) - 1].re = loopFiltOut;
    output->vector.data[emlrtDynamicBoundsCheckR2012b(loop_ub + 1, 1,
      output->size[0], &emlrtBCI, (emlrtConstCTX)sp) - 1].im = DDSOut;
    loopFiltOut = phErr * *integratorGain + *loopFiltState;
    *loopFiltState = loopFiltOut;
    DDSOut = *DDSPreviousInp + *integFiltState;
    *integFiltState = DDSOut;
    *DDSPreviousInp = phErr * *proportionalGain + loopFiltOut;
    *phase = *digitalSynthesizerGain * DDSOut;
    phaseCorrection.vector.data[emlrtDynamicBoundsCheckR2012b(loop_ub + 1, 1,
      phaseCorrection.size[0], &b_emlrtBCI, (emlrtConstCTX)sp) - 1].re = *phase;
    phaseCorrection.vector.data[emlrtDynamicBoundsCheckR2012b(loop_ub + 1, 1,
      phaseCorrection.size[0], &b_emlrtBCI, (emlrtConstCTX)sp) - 1].im = 0.0;
    previousSample->re = output->vector.data[emlrtDynamicBoundsCheckR2012b
      (loop_ub + 1, 1, output->size[0], &c_emlrtBCI, (emlrtConstCTX)sp) - 1].re;
    previousSample->im = output->vector.data[emlrtDynamicBoundsCheckR2012b
      (loop_ub + 1, 1, output->size[0], &c_emlrtBCI, (emlrtConstCTX)sp) - 1].im;
  }

  st.site = &bb_emlrtRSI;
  b = expFunc(*phaseOffset);
  array_creal_T_SetSize(output, output->size[0]);
  loop_ub = output->size[0];
  for (i = 0; i < loop_ub; i++) {
    loopFiltOut = output->vector.data[i].re * b.re - output->vector.data[i].im *
      b.im;
    DDSOut = output->vector.data[i].re * b.im + output->vector.data[i].im * b.re;
    output->vector.data[i].re = loopFiltOut;
    output->vector.data[i].im = DDSOut;
  }

  phaseEstimate_size[0] = phaseCorrection.size[0];
  loop_ub = phaseCorrection.size[0];
  for (i = 0; i < loop_ub; i++) {
    phaseEstimate_data[i] = -(phaseCorrection.vector.data[i].re + *phaseOffset);
  }

  array_creal_T_Destructor(&phaseCorrection);
}

static void array_creal_T_2D_SetSize(coder_array_creal_T_2D *coderArray, int32_T
  size0, int32_T size1)
{
  creal_T *newData;
  int32_T newCapacity;
  int32_T newNumel;
  coderArray->size[0] = size0;
  coderArray->size[1] = size1;
  newNumel = coderArray->size[0] * coderArray->size[1];
  if (newNumel > coderArray->vector.allocated) {
    newCapacity = coderArray->vector.allocated;
    if (newCapacity < 16) {
      newCapacity = 16;
    }

    while (newCapacity < newNumel) {
      if (newCapacity > 1073741823) {
        newCapacity = MAX_int32_T;
      } else {
        newCapacity <<= 1;
      }
    }

    newData = (creal_T *)emlrtMallocMex(sizeof(creal_T) * (uint32_T)newCapacity);
    if (coderArray->vector.data != NULL) {
      memcpy(newData, coderArray->vector.data, sizeof(creal_T) * (uint32_T)
             coderArray->vector.numel);
      if (coderArray->vector.owner) {
        emlrtFreeMex(coderArray->vector.data);
      }
    }

    coderArray->vector.data = newData;
    coderArray->vector.allocated = newCapacity;
    coderArray->vector.owner = true;
  }

  coderArray->vector.numel = newNumel;
}

static void array_creal_T_SetSize(coder_array_creal_T *coderArray, int32_T size0)
{
  creal_T *newData;
  int32_T newCapacity;
  int32_T newNumel;
  coderArray->size[0] = size0;
  newNumel = coderArray->size[0];
  if (newNumel > coderArray->vector.allocated) {
    newCapacity = coderArray->vector.allocated;
    if (newCapacity < 16) {
      newCapacity = 16;
    }

    while (newCapacity < newNumel) {
      if (newCapacity > 1073741823) {
        newCapacity = MAX_int32_T;
      } else {
        newCapacity <<= 1;
      }
    }

    newData = (creal_T *)emlrtMallocMex(sizeof(creal_T) * (uint32_T)newCapacity);
    if (coderArray->vector.data != NULL) {
      memcpy(newData, coderArray->vector.data, sizeof(creal_T) * (uint32_T)
             coderArray->vector.numel);
      if (coderArray->vector.owner) {
        emlrtFreeMex(coderArray->vector.data);
      }
    }

    coderArray->vector.data = newData;
    coderArray->vector.allocated = newCapacity;
    coderArray->vector.owner = true;
  }

  coderArray->vector.numel = newNumel;
}

static void array_creal_T_2D_Constructor(coder_array_creal_T_2D *coderArray)
{
  coderArray->vector.data = (creal_T *)NULL;
  coderArray->vector.numel = 0;
  coderArray->vector.allocated = 0;
  coderArray->vector.owner = true;
  coderArray->size[0] = 0;
  coderArray->size[1] = 0;
}

static void array_creal_T_2D_Destructor(coder_array_creal_T_2D *coderArray)
{
  if (coderArray->vector.owner && (coderArray->vector.data != (creal_T *)NULL))
  {
    emlrtFreeMex(coderArray->vector.data);
  }
}

static void array_creal_T_Constructor(coder_array_creal_T *coderArray)
{
  coderArray->vector.data = (creal_T *)NULL;
  coderArray->vector.numel = 0;
  coderArray->vector.allocated = 0;
  coderArray->vector.owner = true;
  coderArray->size[0] = 0;
}

static void array_creal_T_Destructor(coder_array_creal_T *coderArray)
{
  if (coderArray->vector.owner && (coderArray->vector.data != (creal_T *)NULL))
  {
    emlrtFreeMex(coderArray->vector.data);
  }
}

static void init_simulink_io_address(InstanceStruct_UIzktXqQH1UJHNoLKJcwx
  *moduleInstance)
{
  moduleInstance->emlrtRootTLSGlobal = (void *)cgxertGetEMLRTCtx
    (moduleInstance->S);
  moduleInstance->u0_data = (creal_T (*)[6175])cgxertGetInputPortSignal
    (moduleInstance->S, 0);
  moduleInstance->u0_sizes = (int32_T (*)[2])cgxertGetCurrentInputPortDimensions
    (moduleInstance->S, 0);
  moduleInstance->y0_data = (creal_T (*)[6175])cgxertGetOutputPortSignal
    (moduleInstance->S, 0);
  moduleInstance->y0_sizes = (int32_T (*)[2])
    cgxertGetCurrentOutputPortDimensions(moduleInstance->S, 0);
  moduleInstance->v = (int32_T (*)[2])cgxertGetDWork(moduleInstance->S, 0);
}

/* CGXE Glue Code */
static void mdlOutputs_UIzktXqQH1UJHNoLKJcwx(SimStruct *S, int_T tid)
{
  InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance =
    (InstanceStruct_UIzktXqQH1UJHNoLKJcwx *)cgxertGetRuntimeInstance(S);
  cgxe_mdl_outputs(moduleInstance);
}

static void mdlInitialize_UIzktXqQH1UJHNoLKJcwx(SimStruct *S)
{
  InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance =
    (InstanceStruct_UIzktXqQH1UJHNoLKJcwx *)cgxertGetRuntimeInstance(S);
  cgxe_mdl_initialize(moduleInstance);
}

static void mdlUpdate_UIzktXqQH1UJHNoLKJcwx(SimStruct *S, int_T tid)
{
  InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance =
    (InstanceStruct_UIzktXqQH1UJHNoLKJcwx *)cgxertGetRuntimeInstance(S);
  cgxe_mdl_update(moduleInstance);
}

static void mdlDerivatives_UIzktXqQH1UJHNoLKJcwx(SimStruct *S)
{
  InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance =
    (InstanceStruct_UIzktXqQH1UJHNoLKJcwx *)cgxertGetRuntimeInstance(S);
  cgxe_mdl_derivative(moduleInstance);
}

static mxArray* getSimState_UIzktXqQH1UJHNoLKJcwx(SimStruct *S)
{
  mxArray* mxSS;
  InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance =
    (InstanceStruct_UIzktXqQH1UJHNoLKJcwx *)cgxertGetRuntimeInstance(S);
  mxSS = (mxArray *) cgxe_mdl_get_sim_state(moduleInstance);
  return mxSS;
}

static void setSimState_UIzktXqQH1UJHNoLKJcwx(SimStruct *S, const mxArray *ss)
{
  InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance =
    (InstanceStruct_UIzktXqQH1UJHNoLKJcwx *)cgxertGetRuntimeInstance(S);
  cgxe_mdl_set_sim_state(moduleInstance, emlrtAlias(ss));
}

static void mdlTerminate_UIzktXqQH1UJHNoLKJcwx(SimStruct *S)
{
  InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance =
    (InstanceStruct_UIzktXqQH1UJHNoLKJcwx *)cgxertGetRuntimeInstance(S);
  cgxe_mdl_terminate(moduleInstance);
  free((void *)moduleInstance);
}

static void mdlEnable_UIzktXqQH1UJHNoLKJcwx(SimStruct *S)
{
  InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance =
    (InstanceStruct_UIzktXqQH1UJHNoLKJcwx *)cgxertGetRuntimeInstance(S);
  cgxe_mdl_enable(moduleInstance);
}

static void mdlDisable_UIzktXqQH1UJHNoLKJcwx(SimStruct *S)
{
  InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance =
    (InstanceStruct_UIzktXqQH1UJHNoLKJcwx *)cgxertGetRuntimeInstance(S);
  cgxe_mdl_disable(moduleInstance);
}

static void mdlStart_UIzktXqQH1UJHNoLKJcwx(SimStruct *S)
{
  InstanceStruct_UIzktXqQH1UJHNoLKJcwx *moduleInstance =
    (InstanceStruct_UIzktXqQH1UJHNoLKJcwx *)calloc(1, sizeof
    (InstanceStruct_UIzktXqQH1UJHNoLKJcwx));
  moduleInstance->S = S;
  cgxertSetRuntimeInstance(S, (void *)moduleInstance);
  ssSetmdlOutputs(S, mdlOutputs_UIzktXqQH1UJHNoLKJcwx);
  ssSetmdlInitializeConditions(S, mdlInitialize_UIzktXqQH1UJHNoLKJcwx);
  ssSetmdlUpdate(S, mdlUpdate_UIzktXqQH1UJHNoLKJcwx);
  ssSetmdlDerivatives(S, mdlDerivatives_UIzktXqQH1UJHNoLKJcwx);
  ssSetmdlTerminate(S, mdlTerminate_UIzktXqQH1UJHNoLKJcwx);
  ssSetmdlEnable(S, mdlEnable_UIzktXqQH1UJHNoLKJcwx);
  ssSetmdlDisable(S, mdlDisable_UIzktXqQH1UJHNoLKJcwx);
  cgxe_mdl_start(moduleInstance);

  {
    uint_T options = ssGetOptions(S);
    options |= SS_OPTION_RUNTIME_EXCEPTION_FREE_CODE;
    ssSetOptions(S, options);
  }
}

static void mdlProcessParameters_UIzktXqQH1UJHNoLKJcwx(SimStruct *S)
{
}

void method_dispatcher_UIzktXqQH1UJHNoLKJcwx(SimStruct *S, int_T method, void
  *data)
{
  switch (method) {
   case SS_CALL_MDL_START:
    mdlStart_UIzktXqQH1UJHNoLKJcwx(S);
    break;

   case SS_CALL_MDL_PROCESS_PARAMETERS:
    mdlProcessParameters_UIzktXqQH1UJHNoLKJcwx(S);
    break;

   case SS_CALL_MDL_GET_SIM_STATE:
    *((mxArray**) data) = getSimState_UIzktXqQH1UJHNoLKJcwx(S);
    break;

   case SS_CALL_MDL_SET_SIM_STATE:
    setSimState_UIzktXqQH1UJHNoLKJcwx(S, (const mxArray *) data);
    break;

   default:
    /* Unhandled method */
    /*
       sf_mex_error_message("Stateflow Internal Error:\n"
       "Error calling method dispatcher for module: UIzktXqQH1UJHNoLKJcwx.\n"
       "Can't handle method %d.\n", method);
     */
    break;
  }
}

mxArray *cgxe_UIzktXqQH1UJHNoLKJcwx_BuildInfoUpdate(void)
{
  mxArray * mxBIArgs;
  mxArray * elem_1;
  mxArray * elem_2;
  mxArray * elem_3;
  mxArray * elem_4;
  mxArray * elem_5;
  mxArray * elem_6;
  mxArray * elem_7;
  mxArray * elem_8;
  mxArray * elem_9;
  mxArray * elem_10;
  mxArray * elem_11;
  mxArray * elem_12;
  mxArray * elem_13;
  mxArray * elem_14;
  mxArray * elem_15;
  mxArray * elem_16;
  mxArray * elem_17;
  mxArray * elem_18;
  mxArray * elem_19;
  mxArray * elem_20;
  mxArray * elem_21;
  mxArray * elem_22;
  mxArray * elem_23;
  mxBIArgs = mxCreateCellMatrix(1,3);
  elem_1 = mxCreateCellMatrix(1,6);
  elem_2 = mxCreateCellMatrix(0,0);
  mxSetCell(elem_1,0,elem_2);
  elem_3 = mxCreateCellMatrix(1,4);
  elem_4 = mxCreateString("addIncludeFiles");
  mxSetCell(elem_3,0,elem_4);
  elem_5 = mxCreateCellMatrix(1,3);
  elem_6 = mxCreateString("<string.h>");
  mxSetCell(elem_5,0,elem_6);
  elem_7 = mxCreateString("math.h");
  mxSetCell(elem_5,1,elem_7);
  elem_8 = mxCreateString("signLib.h");
  mxSetCell(elem_5,2,elem_8);
  mxSetCell(elem_3,1,elem_5);
  elem_9 = mxCreateCellMatrix(1,3);
  elem_10 = mxCreateString("");
  mxSetCell(elem_9,0,elem_10);
  elem_11 = mxCreateString("");
  mxSetCell(elem_9,1,elem_11);
  elem_12 = mxCreateString("");
  mxSetCell(elem_9,2,elem_12);
  mxSetCell(elem_3,2,elem_9);
  elem_13 = mxCreateCellMatrix(1,3);
  elem_14 = mxCreateString("");
  mxSetCell(elem_13,0,elem_14);
  elem_15 = mxCreateString("");
  mxSetCell(elem_13,1,elem_15);
  elem_16 = mxCreateString("");
  mxSetCell(elem_13,2,elem_16);
  mxSetCell(elem_3,3,elem_13);
  mxSetCell(elem_1,1,elem_3);
  elem_17 = mxCreateCellMatrix(0,0);
  mxSetCell(elem_1,2,elem_17);
  elem_18 = mxCreateCellMatrix(0,0);
  mxSetCell(elem_1,3,elem_18);
  elem_19 = mxCreateCellMatrix(0,0);
  mxSetCell(elem_1,4,elem_19);
  elem_20 = mxCreateCellMatrix(0,0);
  mxSetCell(elem_1,5,elem_20);
  mxSetCell(mxBIArgs,0,elem_1);
  elem_21 = mxCreateCellMatrix(1,1);
  elem_22 = mxCreateString("comm.internal.utilities.mathCore");
  mxSetCell(elem_21,0,elem_22);
  mxSetCell(mxBIArgs,1,elem_21);
  elem_23 = mxCreateCellMatrix(1,0);
  mxSetCell(mxBIArgs,2,elem_23);
  return mxBIArgs;
}

mxArray *cgxe_UIzktXqQH1UJHNoLKJcwx_fallback_info(void)
{
  const char* fallbackInfoFields[] = { "fallbackType", "incompatiableSymbol" };

  mxArray* fallbackInfoStruct = mxCreateStructMatrix(1, 1, 2, fallbackInfoFields);
  mxArray* fallbackType = mxCreateString("thirdPartyLibs");
  mxArray* incompatibleSymbol = mxCreateString(
    "comm.internal.utilities.mathCore");
  mxSetFieldByNumber(fallbackInfoStruct, 0, 0, fallbackType);
  mxSetFieldByNumber(fallbackInfoStruct, 0, 1, incompatibleSymbol);
  return fallbackInfoStruct;
}
