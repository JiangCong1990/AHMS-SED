module module_rt_rt2d_data

! Qian, to define rt2d for domain scale channel routing
  INTEGER, PARAMETER :: max_domain = 5

#include "rt2d_field_include.inc"

  TYPE(rt2d_field) :: rt2d(max_domain)
  SAVE rt2d

end module module_rt_rt2d_data
