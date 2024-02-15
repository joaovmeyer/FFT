# FFT
Simple implementation of the FFT algorithm and it's uses in polynomial multiplication

A quick overview of the algorithm can be found in FFT.h, maybe I'll add a nicer one here in the future

# Speed Comparison
![345673905734](https://github.com/joaovmeyer/FFT/assets/144701021/325d3b57-2dec-4a29-97cf-ae4621043d80)    
![3945630975634](https://github.com/joaovmeyer/FFT/assets/144701021/bb794dd3-4d0e-46ec-b991-89179f205719)

Here's a execution time comparison of the naive algorithm and the FFT (left image was cut down so we could see the FFT line better). The black line represents the naive algorithm, and the blue the FFT. In the start, the FFT is even a little worse, due to all the overhead it brings, but it shows to be extremelly faster with bigger polynomials, mantaining fast execution times while the naive algorithm skyrockets. At the end (polynomial degree is about 100000), it's clear the FFT can be a powerfull algorithm. My implementation is far from optimal and the FFT proves to be very efficient
