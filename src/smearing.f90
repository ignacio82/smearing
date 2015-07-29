Module smearingmodule
IMPLICIT NONE
contains
subroutine smearing(i, nMCd, nPeriodT, nPeriods, DF, theta, wm) bind(C, name="smearing_")
    use, intrinsic                         :: iso_c_binding, only : c_double, c_int
    integer(c_int), intent(in)                                 :: i, nMCd, nPeriods
    integer(c_int), intent(in), dimension(nPeriods)            :: nPeriodT
    real(c_double), intent(in), dimension(i,5)        :: DF
    real(c_double), intent(in), dimension(i,nMCd)     :: theta
    real(c_double), dimension(nMCd, nPeriods), intent(out)  :: wm
    integer                                             :: d

    do d=1,nMCd
        CALL subsmearing(d, i, nMCd, nPeriodT, nPeriods, DF, theta, wm(d,:))
    end do

end subroutine smearing


subroutine subsmearing(d, i, nMCd, nPeriodT, nPeriods, DF, theta, wm)
    integer, intent(in)                                 :: d, i, nMCd, nPeriods
    integer, intent(in), dimension(nPeriods)            :: nPeriodT
    double precision, intent(in), dimension(i,5)        :: DF
    double precision, intent(in), dimension(i,nMCd)     :: theta
    double precision, dimension(nPeriods), intent(out)  :: wm
    double precision, allocatable, dimension(:,:)       :: out
    double precision, dimension(i,8)                    :: bigDF
    double precision, dimension(i)                      :: epredC, epredT, diff, A1, A2
    double precision, dimension(i,i)                    :: B1, B2
    integer                                             :: j, Period, jj
    double precision                                    :: sumWeights_inv

    diff = 0.0d0
    wm = 0.0d0
    epredC = exp(DF(:,4) + (theta(:,d) * DF(:,1)))
    epredT = exp(DF(:,4) - (theta(:,d) * (1.d0-DF(:,1))))

    do j=1,i
        B1(:,j)=epredT(j)
        B2(:,j)=epredC(j)
    end do

    A1 = matmul(DF(:,5), B1)
    A2 = matmul(DF(:,5), B2)

    diff = (A1-A2)/DBLE(i)

    bigDF(:,1:5) = DF
    bigDF(:,6) = epredC
    bigDF(:,7) = epredT
    bigDF(:,8) = diff

    do jj =1, nPeriods
        ALLOCATE(out(nPeriodT(jj),7))
        CALL subsetPeriod(bigDF, i, 8, nPeriodT(jj), jj, out, sumWeights_inv) !Filters data
        do j=1, nPeriodT(jj)
            wm(jj) = wm(jj) + (out(j,7)*out(j,2)*sumWeights_inv)
        end do
        DEALLOCATE(out)
    end do

end subroutine subsmearing

subroutine subsetPeriod(A, rowA, colA, rowB, Period, B, sumWeights_inv)
    double precision, dimension(rowA, colA), intent(in)    :: A
    integer, intent(in)                                    :: rowA, colA, rowB, Period
    double precision, dimension(rowB,colA-1), intent(out)  :: B
    double precision, intent(out)                          :: sumWeights_inv
    integer                                                :: i, pos
    double precision                                       :: sumWeights

    pos = 1
    do i = 1, size(A,1)
        if(A(i,2)==Period)then
            B(pos,1) = A(i,1)
            B(pos,2:) = A(i,3:)
            pos = pos+1
        end if
    end do
    sumWeights = sum(B(:,2))
    sumWeights_inv = 1.d0/sumWeights
end subroutine subsetPeriod

end module smearingmodule
