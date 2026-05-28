! ***********************************************************************
!
!   Copyright (C) 2010-2025  Bill Paxton & The MESA Team
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU Lesser General Public License
!   as published by the Free Software Foundation,
!   either version 3 of the License, or (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Lesser General Public License for more details.
!
!   You should have received a copy of the GNU Lesser General Public License
!   along with this program. If not, see <https://www.gnu.org/licenses/>.
!
! ***********************************************************************

module run_star_extras

   use star_lib
   use star_def
   use const_def
   use math_lib

   implicit none

   ! these routines are called by the standard run_star check_model
contains

      subroutine extras_controls(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! this is the place to set any procedure pointers you want to change
         ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)
         s% other_neu => tutorial_other_neu

         ! the extras functions in this file will not be called
         ! unless you set their function pointers as done below.
         ! otherwise we use a null_ version which does nothing (except warn).

         s% extras_startup => extras_startup
         s% extras_start_step => extras_start_step
         s% extras_check_model => extras_check_model
         s% extras_finish_step => extras_finish_step
         s% extras_after_evolve => extras_after_evolve
         s% how_many_extra_history_columns => how_many_extra_history_columns
         s% data_for_extra_history_columns => data_for_extra_history_columns
         s% how_many_extra_profile_columns => how_many_extra_profile_columns
         s% data_for_extra_profile_columns => data_for_extra_profile_columns

         s% how_many_extra_history_header_items => how_many_extra_history_header_items
         s% data_for_extra_history_header_items => data_for_extra_history_header_items
         s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
         s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

      end subroutine extras_controls


      subroutine extras_startup(id, restart, ierr)
         integer, intent(in) :: id
         logical, intent(in) :: restart
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_startup


      integer function extras_start_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         extras_start_step = 0
      end function extras_start_step


      ! returns either keep_going, retry, or terminate.
      integer function extras_check_model(id)
     
         use chem_def, only : i_burn_ne, category_name

         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         integer :: i_burn_max

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_check_model = keep_going

         ! if you want to check multiple conditions, it can be useful
         ! to set a different termination code depending on which
         ! condition was triggered.  MESA provides 9 customizable
         ! termination codes, named t_xtra1 .. t_xtra9.  You can
         ! customize the messages that will be printed upon exit by
         ! setting the corresponding termination_code_str value.
         ! termination_code_str(t_xtra1) = 'my termination condition'

         ! determine the category of maximum burning
         i_burn_max = maxloc(s% L_by_category,1)

         ! stop if the luminosity is dominated by neon burning
         if ( i_burn_max .eq. i_burn_ne) then
            extras_check_model = terminate
            s% termination_code = t_xtra1
            termination_code_str(t_xtra1) = 'neon burning is dominant'
            return
         end if

         ! by default, indicate where (in the code) MESA terminated
         if (extras_check_model == terminate) s% termination_code = t_extras_check_model
        

      end function extras_check_model


      integer function how_many_extra_history_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_columns = 2
      end function how_many_extra_history_columns


      subroutine data_for_extra_history_columns(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s

         real(dp), parameter :: frac = 0.90
         integer :: i
         real(dp) :: edot, edot_partial

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! calculate the total nuclear energy release rate by integrating
         ! the specific rate (eps_nuc) over the star.  using the dot_product
         ! intrinsic is a common idiom for calculating integrated quantities.
         ! note that one needs to explicitly limit the range of the arrays.
         ! NEVER assume that the array size is the same as the number of zones.
         edot = dot_product(s% dm(1:s% nz), s% eps_nuc(1:s% nz))

         ! the center of the star is at i = s% nz and the surface at i = 1 .
         ! so go from the center outward until 90% of the integrated eps_nuc
         ! is enclosed.  exit and then i will contain the desired cell index.
         edot_partial = 0
         do i = s% nz, 1, -1
            edot_partial = edot_partial + s% dm(i) * s% eps_nuc(i)
            if (edot_partial .ge. (frac * edot)) exit
         end do

         ! note: do NOT add these names to history_columns.list
         ! the history_columns.list is only for the built-in log column options.
         ! it must not include the new column names you are adding here.

         ! column 1
         names(1) = "m90"
         vals(1) = s% q(i) * s% star_mass  ! in solar masses

         ! column 2
         names(2) = "log_R90"
         vals(2) = log10(s% R(i) / rsun) ! in solar radii

         ierr = 0
      end subroutine data_for_extra_history_columns


      integer function how_many_extra_profile_columns(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_columns = 1 
      end function how_many_extra_profile_columns


      subroutine data_for_extra_profile_columns(id, n, nz, names, vals, ierr)
         integer, intent(in) :: id, n, nz
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(nz,n)
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         integer :: k
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         ! note: do NOT add the extra names to profile_columns.list
         ! the profile_columns.list is only for the built-in profile column options.
         ! it must not include the new column names you are adding here.

         ! here is an example for adding a profile column
         if (n /= 1) stop 'data_for_extra_profile_columns'
         names(1) = 'beta'
         do k = 1, nz
            vals(k,1) = s% Pgas(k)/s% Peos(k)
         end do

      end subroutine data_for_extra_profile_columns


      integer function how_many_extra_history_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_history_header_items = 0
      end function how_many_extra_history_header_items


      subroutine data_for_extra_history_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_history_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra history header item
         ! also set how_many_extra_history_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_history_header_items


      integer function how_many_extra_profile_header_items(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
         how_many_extra_profile_header_items = 0
      end function how_many_extra_profile_header_items


      subroutine data_for_extra_profile_header_items(id, n, names, vals, ierr)
         integer, intent(in) :: id, n
         character (len=maxlen_profile_column_name) :: names(n)
         real(dp) :: vals(n)
         type(star_info), pointer :: s
         integer, intent(out) :: ierr
         ierr = 0
         call star_ptr(id,s,ierr)
         if(ierr/=0) return

         ! here is an example for adding an extra profile header item
         ! also set how_many_extra_profile_header_items
         ! names(1) = 'mixing_length_alpha'
         ! vals(1) = s% mixing_length_alpha

      end subroutine data_for_extra_profile_header_items


      ! returns either keep_going or terminate.
      ! note: cannot request retry; extras_check_model can do that.
      integer function extras_finish_step(id)
         integer, intent(in) :: id
         integer :: ierr
         type (star_info), pointer :: s

         integer :: f

         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return

         extras_finish_step = keep_going

         ! MESA provides a number of variables that make it easy to get user input.
         ! these are part of the star_info structure and are named
         ! x_character_ctrl, x_integer_ctrl, x_logical_ctrl, and x_ctrl.
         ! by default there are num_x_ctrls, which defaults to 100, of each.
         ! they can be specified in the controls section of your inlist.

         f = s% x_integer_ctrl(1)

         ! MESA also provides a number of arrays that are useful for implementing
         ! algorithms which require a state. if you use these variables
         ! restarts and retries will work without doing anything special.
         ! they are named xtra, ixtra, lxtra.
         ! they are automatically versioned, that is if you set s% xtra(1), then
         ! s% xtra_old(1) will contains the value of s% xtra(1) from the previous step.

         s% xtra(1) = s% log_center_density

         ! this expression will evaluate to true if f times the log center density
         ! has crossed an integer during the last step.  If f = 5, then we will get
         ! output at log center density = {... 1.0, 1.2, 1.4, 1.6, 1.8, 2.0 ... }
         if ((floor(f * s% xtra_old(1)) - floor(f * s% xtra(1)) .ne. 0)) then

            ! save a profile & update the history
            s% need_to_update_history_now = .true.
            s% need_to_save_profiles_now = .true.

            ! by default the priority is 1; you can change that if you'd like
            ! s% save_profiles_model_priority = ?

         endif

         ! see extras_check_model for information about custom termination codes
         ! by default, indicate where (in the code) MESA terminated
         if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

      end function extras_finish_step


      subroutine extras_after_evolve(id, ierr)
         integer, intent(in) :: id
         integer, intent(out) :: ierr
         type (star_info), pointer :: s
         ierr = 0
         call star_ptr(id, s, ierr)
         if (ierr /= 0) return
      end subroutine extras_after_evolve

      subroutine tutorial_other_neu(  &
           id, k, T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim, flags, &
           loss, sources, ierr)
         use neu_lib, only: neu_get
         use neu_def
         integer, intent(in) :: id ! id for star
         integer, intent(in) :: k ! cell number or 0 if not for a particular cell
         real(dp), intent(in) :: T ! temperature
         real(dp), intent(in) :: log10_T ! log10 of temperature
         real(dp), intent(in) :: Rho ! density
         real(dp), intent(in) :: log10_Rho ! log10 of density
         real(dp), intent(in) :: abar ! mean atomic weight
         real(dp), intent(in) :: zbar ! mean charge
         real(dp), intent(in) :: log10_Tlim
         logical, intent(inout) :: flags(num_neu_types) ! true if should include the type of loss
         real(dp), intent(inout) :: loss(num_neu_rvs) ! total from all sources
         real(dp), intent(inout) :: sources(num_neu_types, num_neu_rvs)
         integer, intent(out) :: ierr

         ! before we can use controls associated with the star we need to get access
         type (star_info), pointer :: s
         call star_ptr(id, s, ierr)
         if (ierr /= 0) then ! OOPS
            return
         end if

         ! separately control whether each type of neutrino loss is included
         flags(pair_neu_type) = s% x_logical_ctrl(1)
         flags(plas_neu_type) = s% x_logical_ctrl(2)
         flags(phot_neu_type) = s% x_logical_ctrl(3)
         flags(brem_neu_type) = s% x_logical_ctrl(4)
         flags(reco_neu_type) = s% x_logical_ctrl(5)

         ! the is the normal routine that MESA provides
         call neu_get(  &
             T, log10_T, Rho, log10_Rho, abar, zbar, log10_Tlim, flags, &
             loss, sources, ierr)

      end subroutine tutorial_other_neu


end module run_star_extras
