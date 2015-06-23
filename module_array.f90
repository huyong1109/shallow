module module_array

    use module_para

    implicit none

    real*8,  dimension(1:n1,1:n) ::    u               ! zonal wind
    real*8,  dimension(1:n1,1:n) ::    v               ! meridional wind
    real*8,  dimension(1:n1,1:n) ::    wh              ! geopotential height

    real*8,  dimension(1:n1,1:n) ::    wu              ! u*sqrt(wh)
    real*8,  dimension(1:n1,1:n) ::    wv              ! v*sqrt(wh)(

end module module_array
!
