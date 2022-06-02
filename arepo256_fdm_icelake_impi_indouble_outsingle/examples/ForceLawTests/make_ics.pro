
N = 100

seed = 42L

pos = randomu(seed, 3, N)

pos *= 1.0
pos += 49.5

mass = fltarr(N)

id = lindgen(N)+1

mass(0) = 1 

vel = fltarr(3, N)

u = fltarr(N)
u(*) = 1



openw,1,"ICs/forcelawtest_ics.dat",/f77_unformatted

npart=lonarr(6)	
massarr=dblarr(6)
time=0.0D
redshift=0.0D
flag_sfr=0L
flag_feedback=0L

bytesleft=256-6*4 - 6*8 - 8 - 8 - 2*4
la=intarr(bytesleft/2)

npart(2)=N

writeu,1,npart,massarr,time,redshift,flag_sfr,flag_feedback,la
writeu,1,pos
writeu,1,vel
writeu,1,id
writeu,1,mass
writeu,1,u
writeu,1,u
close,1


ende:
end

