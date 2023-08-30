module SPT3G_utils

  implicit None
#ifdef _STANDALONE_

  ! minimal code, lifted from CAMB et al to make the thing standalone
  
  integer, parameter :: mcp= KIND(1.d0)
  double precision, parameter :: pi=3.14159265358979323846264338328d0

  Type :: TStringList
    integer :: Count =0
    integer :: max =0
    character(len = 1024),dimension(:),allocatable :: Item
    contains
      procedure :: SetFromString
      procedure :: add      
  end type TStringList

  contains
    subroutine add(this,string1024)
      class(TStringList) :: this
      integer :: oldmax,i
      character(len=*), intent(in)::string1024
      character(len=1024),dimension(:),allocatable :: oldItem
      if (this%max==0) then
        this%max=16
        allocate(this%Item((this%max)))      
      endif
      if (this%max==this%Count) then
        allocate(oldItem(this%max))
        do i=1,this%max
          oldItem(i) = trim(this%Item(i))
        enddo
        deallocate(this%Item)

        allocate(this%Item((this%max)*2))
        do i=1,this%max
          this%Item(i) = trim(oldItem(i))
        enddo
        deallocate(oldItem)
        this%max = this%max*2
      endif
      this%Count = this%Count + 1
      !print *,"'",string1024,"'",len(string1024),size(this%Item),len(this%Item(1)),this%Count,this%max
      this%Item(this%Count) = trim(string1024(:len(string1024)))
      !print *,"it:'",this%Item(this%Count),"'","'",trim(this%Item(this%Count)),"'"
    end subroutine

    subroutine SetFromString(this, longstring)
      class(TStringList) :: this
      character(Len=*), intent(in) :: longstring
      character(LEN=1024), allocatable :: item
      integer i,j
      character(LEN=:), allocatable :: valid_chars

      valid_chars='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-.'
      allocate(item)
      j=0
      do i=1, len_trim(longstring)
          if (verify(longstring(i:i),valid_chars) == 0) then
              j=j+1
              item(j:j) = longstring(i:i)
          else
              if (trim(longstring(i:i))/='') then
                  write (*,*) 'Invalid character in: '//trim(longstring)
              end if
              if (j>0) call this%Add(item(1:j))
              j=0
          end if
      end do
      if (j>0) call this%Add(item(1:j))

      
    end subroutine SetFromString
#endif

end module
