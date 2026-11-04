!> Separate compilation unit avoids a backend/helper/MPI module cycle.
subroutine s_poem_abort(message)

    use m_mpi_proxy, only: s_mpi_abort

    implicit none

    character(len=*), intent(in) :: message

    call s_mpi_abort(message)

end subroutine s_poem_abort
