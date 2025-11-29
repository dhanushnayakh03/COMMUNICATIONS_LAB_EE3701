/* Include files */

#include "plutoradioQPSKReceiverSimulinkExample_cgxe.h"
#include "m_Vy8xGh9WJGSVoTt911H7G.h"

unsigned int cgxe_plutoradioQPSKReceiverSimulinkExample_method_dispatcher
  (SimStruct* S, int_T method, void* data)
{
  if (ssGetChecksum0(S) == 2193679718 &&
      ssGetChecksum1(S) == 188094155 &&
      ssGetChecksum2(S) == 2403564237 &&
      ssGetChecksum3(S) == 1747372316) {
    method_dispatcher_Vy8xGh9WJGSVoTt911H7G(S, method, data);
    return 1;
  }

  return 0;
}
