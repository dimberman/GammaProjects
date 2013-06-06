#include <stdio.h>
#include "Gamma/AudioIO.h"
#include "Gamma/Oscillator.h"
#include "../examples.h"
#include <math.h>
#include "assert.h"


using namespace std;
//using namespace al;
//using namespace osc;
float bpmInput=20;




//a reverb class ripped from allocore. 

template <class T>
class Reverb{
public:
    
	Reverb()
	:	mPreDelay(10),
    mAPIn1(142), mAPIn2(107), mAPIn3(379), mAPIn4(277),
    mAPDecay11(672), mAPDecay12(1800), mDly11(4453), mDly12(3720),
    mAPDecay21(908), mAPDecay22(2656), mDly21(4217), mDly22(3163)
	{
        //		bandwidth(0.9995);
        //		decay(0.5);
        //		damping(0.0005);
        //		diffusion(0.75, 0.625, 0.7, 0.5);
		bandwidth(0.9995);
		decay(0.85);
		damping(0.4);
		diffusion(0.76, 0.666, 0.707, 0.571);
	}
    
	/// Set input signal bandwidth, [0,1]
	Reverb& bandwidth(T v){ mOPIn.damping(T(1)-v); return *this; }
    
	/// Set high-frequency damping amount, [0, 1]
	Reverb& damping(T v){ mOP1.damping(v); mOP2.damping(v); return *this; }
    
	/// Set decay rate, [0, 1)
	Reverb& decay(T v){ mDecay=v; return *this; }
    
	/// Set diffusion amounts
	
	/// The recommended range of these coefficients is from 0.0 to 0.9999999
	///
	Reverb& diffusion(T in1, T in2, T decay1, T decay2){
		mDfIn1=in1;	mDfIn2=in2; mDfDcy1=decay1; mDfDcy2=decay2;
		return *this;
	}
	
	/// Set input diffusion 1 amount, [0,1)
	Reverb& diffusionIn1(const T& v){ mDfIn1=v; return *this; }
	
	/// Set input diffusion 2 amount, [0,1)
	Reverb& diffusionIn2(const T& v){ mDfIn2=v; return *this; }
	
	/// Set tank decay diffusion 1 amount, [0,1)
	Reverb& diffusionDecay1(const T& v){ mDfDcy1=v; return *this; }
	
	/// Set tank decay diffusion 2 amount, [0,1)
	Reverb& diffusionDecay2(const T& v){ mDfDcy2=v; return *this; }
    
	/// Compute wet stereo output from dry mono input
	void operator()(const T& i0, T& o1, T& o2, T gain = T(0.6)){
		T v = mPreDelay(i0 * T(0.5));
		v = mOPIn(v);
		v = mAPIn1.comb(v, mDfIn1,-mDfIn1);
		v = mAPIn2.comb(v, mDfIn1,-mDfIn1);
		v = mAPIn3.comb(v, mDfIn2,-mDfIn2);
		v = mAPIn4.comb(v, mDfIn2,-mDfIn2);
		
		T a = v + mDly22.back() * mDecay;
		T b = v + mDly12.back() * mDecay;
		
		a = mAPDecay11.comb(a,-mDfDcy1, mDfDcy1);
		a = mDly11(a);
		a = mOP1(a) * mDecay;
		a = mAPDecay12.comb(a, mDfDcy2,-mDfDcy2);
		mDly12.write(a);
        
		b = mAPDecay21.comb(b,-mDfDcy1, mDfDcy1);
		b = mDly21(b);
		b = mOP2(b) * mDecay;
		b = mAPDecay22.comb(b, mDfDcy2,-mDfDcy2);
		mDly22.write(b);
		
		o1 = gain*(  mDly21.read(266)
                   + mDly21.read(2974)
                   - mAPDecay22.read(1913)
                   + mDly22.read(1996)
                   - mDly11.read(1990)
                   - mAPDecay12.read(187)
                   - mDly12.read(1066));
        
		o2 = gain*(  mDly11.read(353)
                   + mDly11.read(3627)
                   - mAPDecay12.read(1228)
                   + mDly12.read(2673)
                   - mDly21.read(2111)
                   - mAPDecay22.read(335)
                   - mDly22.read(121));
	}
    
	/// Compute wet/dry mix stereo output from dry mono input
	
	/// \returns dry sample
	///
	T mix(T& io0, T& o1, T wetAmt){
		T s = io0;
		(*this)(s, io0, o1, wetAmt*T(0.6));
		io0 += s;
		o1  += s;
		return s;
	}
    
protected:
	class DelayLine {
	public:
		DelayLine(int size)
		:	mPos(0), mSize(0), mBuf(0)
		{	resize(size); }
        
		~DelayLine(){ deleteBuf(); }
        
		/// Read value at delay i
		const T& read(int i) const {
			int ind = pos()-i;
			if(ind < 0) ind += size();
			//else if(ind >= size()) ind -= size();
			return mBuf[ind];
		}
		
		/// Write value to delay
		void write(const T& v){
			mBuf[pos()] = v;
			++mPos; if(mPos >= size()) mPos=0;
		}
        
        
		const T& back() const { return mBuf[indexBack()]; }
        
		int indexBack() const {
			int i = pos()+1;
			return (i < size()) ? i : 0;
		}
        
		/// Get absolute index of write tap
		int pos() const { return mPos; }
		
		int size() const { return mSize; }
        
        
		
		/// Write new value and return oldest value
		T operator()(const T& v){
			T r = mBuf[pos()];
			write(v);
			return r;
		}
		
		T comb(const T& v, const T& ffd, const T& fbk){
			T d = mBuf[pos()];
			T r = v + d*fbk;
			write(r);
			return d + r*ffd;
		}
        
		void resize(int n){
			if(n != mSize){
				deleteBuf();
				mBuf = (T*)::calloc(n, sizeof(T));
				mSize = n;
				if(mPos >= n) mPos = mPos % n;
			}
		}
        
	protected:
		void deleteBuf(){ if(mBuf) ::free(mBuf); mBuf=0; }
        
		int mPos;
		int mSize;
		T * mBuf;
	};
    
	class OnePole{
	public:
		OnePole(): mO1(0), mA0(1), mB1(0){}
		void damping(const T& v){ mB1=v; mA0=T(1)-v; }
		T operator()(const T& i0){ return mO1 = mO1*mB1 + i0*mA0; }
	protected:
		T mO1;
		T mA0, mB1;
	};
    
	T mDfIn1, mDfIn2, mDfDcy1, mDfDcy2, mDecay;
	DelayLine mPreDelay;
	OnePole mOPIn;
	DelayLine mAPIn1, mAPIn2, mAPIn3, mAPIn4;
	DelayLine mAPDecay11, mAPDecay12, mDly11, mDly12;
	OnePole mOP1;
	DelayLine mAPDecay21, mAPDecay22, mDly21, mDly22;
	OnePole mOP2;
};









































#include <stdio.h>
#include "Gamma/AudioIO.h"
#include "Gamma/Oscillator.h"
#include "../examples.h"

//#include <Noise.h>
//#include <effects.h>
using namespace gam;


float beat=bpmInput/60;
//float bpm=60;
//float beat = bpm/60;

//Chirp<> src(40, 20, 0.01);
SineD<> src2(440,1,.3,0);
SineD<> src3(420,1,.5,0);
SineD<> src4(610,1,.5,0);
SineD<> src5(260,1,.5,0);
SineD<> src6(153,1,.5,0);
SineD<> src(802,1, .5,  0);
//first beat

Accum<> tmr2(beat,.335*beat);
Accum<> tmr6(beat,.137*beat);
Accum<> tmr7(beat,.240*beat);

//second beat
Accum<> tmr(beat);
Accum<> tmr3(beat,.35*beat);
Accum<> tmr4(beat,.101*beat);
Accum<> tmr5(beat,.602*beat);
Accum<> sweep(10, 0);
Accum<> tmrRand(.05, 0);


//Comb<> cmb1(1, -.9, .9);
AD<> env(0.01, .03, 1);
AD<> env2(0.01, .03, 1);
AD<> envLPF(0.01, .03, 1);
AD<> envwhite(1,2, 1);
AD<> envswell(4,4,2);
Curve<> envpan(434200,-3);
Curve<> highDropPan(5,-3,1,.5);
Biquad<> LPF(120, 20, LOW_PASS);
Biquad<> BPF(130, 20, BAND_PASS);
NoiseWhite<>swnoise(2);
//Biquad<> HPF(30, 20, BAND_PASS);
NoisePink<> noise(3);
#define FD 1
Comb<> cmb1(1, -FD, FD), cmb2(1, -FD, FD), cmb3(1, -FD, FD);

//Blob b;

int port = 11111;
const char * addr = "localhost";


struct position {
    double x, y, z;
};
float hpfreq=80.0;

Reverb<float> reverb;
void audioCB(AudioIOData& io){
		position pos;
		pos.x=-10;
		pos.y=0;
		pos.z=0;
	while(io()){
		float r;
		r=rand()%7;

	//	cout<<bpmInput<<endl;
	//	cout<<tmr2.freq()<<endl;
		float s = 0;
		float tone = 0;
		float bass = 0;
		if (tmr()){	

			src3.reset();
			envLPF.reset();
			}
	if(tmr3()){
			env.reset();
			env2.reset();
			//src3.reset();
			tmr3.freq(rnd::uni(beat,beat*2));
			src2.reset();
		}
		if(tmr2()){
//			src.reset();
			tmr2.freq(rnd::uni(beat,beat*3));
			src3.reset();
			src2.reset();
		}
		if (tmr4()){	
	//		src.reset();
			pos.x+=1;
			tmr4.freq(rnd::uni(beat,beat*2));
			src3.reset();
			envLPF.reset();
			}
		if (tmr5()){	
	//		src.reset();
			pos.x+=1;
			tmr5.freq(rnd::uni(beat/2,beat*2));
			src4.reset();
			envLPF.reset();
			}
		if (tmr6()){	
	//		src.reset();
			pos.x+=1;
			src5.reset();
			envLPF.reset();
			}
		if (tmr7()){	
	//		src.reset();
			pos.x+=1;
			src6.reset();
			envLPF.reset();
			}
		if (tmrRand())
		{
			envswell.reset();
			envpan.set(434200,-3);
//			beat*=(float)r;

		}
		if(sweep()){

				envwhite.reset();
		

		}
			
		float wet1, wet2;
		s=src()*env()*.3;
		hpfreq*=envLPF();
		tone = src2()*.2;
	 	bass = src3()*.3;
	 	s+=src4()*.2+src5()*.2+src6()*.2;
		tone+=bass;
		s+=tone;
		reverb(s,wet1, wet2);
		//double xplane =  pos.x/10;




		//define my sound variables with pan
		float batleft=swnoise()*envwhite()*envswell()*.15*envpan();
		float batright= swnoise()*envwhite()*(envswell())*.15*(1-envpan());









		float preout=LPF(s);
		float out = preout+LPF(wet1*.2);
		io.out(0) = out+batleft;
		io.out(1) = out+batright;
	}
}

int main (int argc, char * argv[]){
        
	reverb.bandwidth(0.9);		// low-pass amount on input
	reverb.damping(0.5);		// high-frequency damping
	reverb.decay(0.8);			// tail decay factor
	reverb.diffusion(0.76, 0.666, 0.707, 0.571); // diffusion amounts

	AudioIO audioIO(256, 44100, audioCB, 0, 2, 1);
		Sync::master().spu(audioIO.fps());
	audioIO.start();
	
	
	//this is me attempting to recieve an osc message
/*
	osc::Send s(port, addr);
	osc::Recv r(port);
	r.handler(handler);
	r.timeout(0.1); // set receiver to block with timeout
	r.start();
	s.beginBundle(1);
	s.beginMessage("/test");
	s<<"notastring";
	s.endMessage();
	s.endBundle();
	s.send();


*/	
	printf("\nPress 'enter' to quit...\n"); getchar();
	return 0;
}