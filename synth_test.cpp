/**
 * @file 	synth_test.cpp
 *
 * @date 	20.11.2022
 * @author 	Martin Helmich <helmima1@fel.cvut.cz>
 * @author 	Veronika Ovsyannikova <ovsyaver@fel.cvut.cz>
 *
 * The coolest program on earth.
 */


#include "iimavlib.h"
#include "iimavlib/WaveSource.h"
#include "iimavlib/WaveSink.h"
#include "iimavlib/Utils.h"
#include "iimavlib/AudioFilter.h"
#include "iimavlib/filters/SineMultiply.h"
#include "iimavlib/filters/SimpleEchoFilter.h"
#include <string>
#include <cassert>
#include <cmath>
#include "iimavlib/SDLDevice.h"
#include "iimavlib_high_api.h"
#include "iimavlib/video_ops.h"
#ifdef SYSTEM_LINUX
#include <unistd.h>
#endif
#include <algorithm>
#include <random>

using namespace iimavlib;

namespace {
#ifdef MODERN_COMPILER
	// mapping of keys to notes (Notes C4 - C6)
	const std::map<char, double> notes = {
	{'a', 261.63}, // C4
	{'s', 293.66}, // D4
	{'d', 329.63}, // E4
	{'f', 349.23}, // F4
	{'g', 392.00}, // G4
	{'h', 440.00}, // A4
	{'j', 493.88}, // B4
	{'k', 523.25}, // C5
	{'l', 587.33}, // D5
	{'z', 659.25}, // E5
	{'x', 698.46}, // F5
	{'c', 783.99}, // G5
	{'v', 880.00}, // A5
	{'b', 987.77}, // B5
	{'n', 1046.50}, // C6
	{'m', 1174.66}, //D6
	};
#else
	const std::map<char, double> notes = iimavlib::InitMap<char, double>
		('a', 261.63) // C4
		('s', 293.66) // D4
		('d', 329.63) // E4
		('f', 349.23) // F4
		('g', 392.00) // G4
		('h', 440.00) // A4
		('j', 493.88) // B4
		('k', 523.25) // C5
		('l', 587.33) // D5
		('z', 659.25) // E5
		('x', 698.46) // F5
		('c', 783.99) // G5
		('v', 880.00) // A5
		('b', 987.77) // B5
		('n', 1046.50) // C6
		('m', 1174.66); //D6
	#endif

namespace {
	// Max value for int16_t
	const double max_val = std::numeric_limits<int16_t>::max();

	// Value of 2*PI
	const double pi2 = 8.0 * std::atan(1.0);
}

}

// ------ GENERATORS ------ 
class NoiseGenerator : public AudioFilter
{
public:
	NoiseGenerator(double frequency) :AudioFilter(pAudioFilter())
	{
	}
	
private:
	error_type_t do_process(audio_buffer_t& buffer)
	{
		std::random_device rd; // obtain a random number from hardware
		std::mt19937 gen(rd()); // seed the generator
		std::uniform_int_distribution<> distr(-10000, 10000); // define the range	
		for (auto& sample : buffer.data) {
			sample = distr(gen);
		}
		return error_type_t::ok;
	}
	
	std::mutex frequency_mutex_;
};

class AddNoise : public AudioFilter
{
public:
	AddNoise() :AudioFilter(pAudioFilter())
	{
	}

private:
	error_type_t do_process(audio_buffer_t& buffer)
	{
		std::random_device rd; // obtain a random number from hardware
		std::mt19937 gen(rd()); // seed the generator
		std::uniform_int_distribution<> distr(-10000, 10000); // define the range	
		for (auto& sample : buffer.data) {
			sample = sample * 0.5 + distr(gen)* 0.5;
		}
		return error_type_t::ok;
	}

	std::mutex frequency_mutex_;
};

enum WaveShape { Sine, Triangle, Saw, Square };

class WaveShapesGeneratorOld : public AudioFilter
{
public:
	WaveShapesGeneratorOld(double frequency, WaveShape shape) : AudioFilter(pAudioFilter()), shape(shape)
	{
	}

	void change_shape(WaveShape shape) {
		this->shape = shape;
	}

	void set_frequency(double frequency) {
		std::unique_lock<std::mutex> lock(frequency_mutex_);

		frequency_ = frequency;
	}

private:
	error_type_t do_process(audio_buffer_t& buffer)
	{
		std::unique_lock<std::mutex> lock(frequency_mutex_);
		const double step = 1.0 / convert_rate_to_int(buffer.params.rate);
		for (auto& sample : buffer.data) {
			sample = static_cast<int16_t>(max_val * std::sin(time_ * frequency_ * pi2));
			time_ = time_ + step;
		}
		buffer.valid_samples = buffer.data.size();
		return error_type_t::ok;
		switch (shape) {
		case WaveShape::Sine:
			
			break;
		/*
		case WaveShape::Triangle:
			break;
		case WaveShape::Saw:
			break;
		case WaveShape::Square:
			break;
			*/
		default:
			break;
		}
	}

	WaveShape shape;
	double frequency_;
	double time_;
	std::mutex frequency_mutex_;
};

/**
* Generator abstract class for Sine, Triangle, Saw and Square wave shapes. In addition to AudioFilter class it has set_frequency method, otherweise it is the same as 
* AudioFilter.
* Reason of making this abstract class is better handling in control.
*/
class WaveGenerator : public AudioFilter
{
public:
	WaveGenerator(const pAudioFilter& child, double frequency_) : AudioFilter(pAudioFilter()), child_(child), frequency_(frequency_), time_(0.0) {}
	/**
	 * @brief Destructor
	 */
	virtual ~WaveGenerator() {}

	void set_frequency(double frequency) {
		std::unique_lock<std::mutex> lock(frequency_mutex_);

		frequency_ = frequency;
	}

	error_type_t process(audio_buffer_t& buffer)
	{
		error_type_t ret = error_type_t::ok;

		if (child_) {
			ret = child_->process(buffer);
			if (ret != error_type_t::ok) return ret;
		}

		return do_process(buffer);
	}

	audio_params_t get_params() const {
		return do_get_params();
	}

	pAudioFilter get_child(size_t depth = 0)
	{
		if (!depth) return child_;
		if (!child_) return pAudioFilter();
		return child_->get_child(depth - 1);
	}
private:

	virtual error_type_t do_process(audio_buffer_t& buffer) = 0;

	virtual audio_params_t do_get_params() const 
	{
		if (child_) {
			return child_->get_params();
		}
		return audio_params_t();
	}
protected: 
	pAudioFilter child_;
	double frequency_;
	double time_;
	std::mutex frequency_mutex_;
};

class WaveShapesGenerator : public AudioFilter
{
public:
	WaveShapesGenerator(double frequency, WaveShape shape) : AudioFilter(pAudioFilter()), shape(shape), frequency_(frequency),
		time_(0.0), sample_prep(0), direction(true), adding_coeff(0)
	{
	}
	void set_frequency(double frequency) {
		std::unique_lock<std::mutex> lock(frequency_mutex_);

		frequency_ = frequency;
	}
private:
	error_type_t do_process(audio_buffer_t& buffer)
	{
		const audio_params_t& params = buffer.params;
		const double step = 1.0 / convert_rate_to_int(buffer.params.rate);
		std::unique_lock<std::mutex> lock(frequency_mutex_);
		switch (shape) {
		case WaveShape::Sine:
			for (auto& sample : buffer.data) {
				sample = static_cast<int16_t>(max_val * std::sin(time_ * frequency_ * pi2));
				time_ = time_ + step;
			}
			break;
		case WaveShape::Triangle:
			adding_coeff = 65536 / (convert_rate_to_int(params.rate) / frequency_);
			for (auto& sample : buffer.data) {
				sample = sample_prep;
				int16_t temp;
				if (direction) {
					temp = sample_prep + adding_coeff;
					if (sample_prep > temp) {
						direction = false;
						sample_prep -= adding_coeff;
					}
					else {
						sample_prep += adding_coeff;
					}
				}
				else {
					temp = sample_prep - adding_coeff;
					if (sample_prep < temp) {
						direction = true;
						sample_prep += adding_coeff;
					}
					else {
						sample_prep -= adding_coeff;
					}
				}
			}
			break;
		case WaveShape::Saw:
			adding_coeff = 65536 / (convert_rate_to_int(params.rate) / frequency_);
			for (auto& sample : buffer.data) {
				sample = sample_prep;
				sample_prep += adding_coeff;
			}
			break;
		case WaveShape::Square:
			for (auto& sample : buffer.data) {
				sample_prep = static_cast<int16_t>(max_val * std::sin(time_ * frequency_ * pi2));
				if (sample_prep > 0) {
					sample = 10000;
				}
				else {
					sample = -10000;
				}
				time_ = time_ + step;
			}
			break;
		default:
			break;
		}
		buffer.valid_samples = buffer.data.size();
		return error_type_t::ok;
	}

	WaveShape shape;
	double frequency_;
	double time_;
	// Triangle wave parameters
	int16_t adding_coeff;
	int16_t sample_prep;
	boolean direction;

	std::mutex frequency_mutex_;
};

class SineGenerator : public WaveGenerator
{
public:
	SineGenerator(double frequency) : WaveGenerator(pAudioFilter(), frequency)
	{
	}
	void set_frequency(double frequency) {
		std::unique_lock<std::mutex> lock(frequency_mutex_);
		
		frequency_ = frequency;
	}
private:
	error_type_t do_process(audio_buffer_t& buffer)
	{
		const double step = 1.0 / convert_rate_to_int(buffer.params.rate);
		std::unique_lock<std::mutex> lock(frequency_mutex_);
		for (auto& sample : buffer.data) {
			sample = static_cast<int16_t>(max_val * std::sin(time_ * frequency_ * pi2));
			time_ = time_ + step;
		}
		buffer.valid_samples = buffer.data.size();
		return error_type_t::ok;
	}

	//double frequency_;
	//double time_;
	std::mutex frequency_mutex_;
};

class SquareWaveGenerator : public WaveGenerator
{
public:
	SquareWaveGenerator(double frequency) :WaveGenerator(pAudioFilter(), frequency)
	{
	}
	/*
	void set_frequency(double frequency) {
		std::unique_lock<std::mutex> lock(frequency_mutex_);

		frequency_ = frequency;
	}
	*/
private:
	error_type_t do_process(audio_buffer_t& buffer)
	{
		const double step = 1.0 / convert_rate_to_int(buffer.params.rate);
		std::unique_lock<std::mutex> lock(frequency_mutex_);
		int16_t sample_prep;
		for (auto& sample : buffer.data) {
			sample_prep = static_cast<int16_t>(max_val * std::sin(time_ * frequency_ * pi2));
			if (sample_prep > 0) {
				sample = 10000;
			}
			else {
				sample = -10000;
			}
			time_ = time_ + step;
		}
		buffer.valid_samples = buffer.data.size();
		return error_type_t::ok;
	}

	//double frequency_;
	//double time_;
	std::mutex frequency_mutex_;
};

class SawWaveGenerator : public WaveGenerator
{
public:
	SawWaveGenerator(double frequency) :WaveGenerator(pAudioFilter(), frequency)
	{
	}
	/*
	void set_frequency(double frequency) {
		std::unique_lock<std::mutex> lock(frequency_mutex_);

		frequency_ = frequency;
	*/
private:
	error_type_t do_process(audio_buffer_t& buffer)
	{
		std::unique_lock<std::mutex> lock(frequency_mutex_);
		const audio_params_t& params = buffer.params;
		adding_coeff = 65536 / (convert_rate_to_int(params.rate) / frequency_);
		for (auto& sample : buffer.data) {
			sample = sample_prep;
			sample_prep += adding_coeff;
		}
		buffer.valid_samples = buffer.data.size();
		return error_type_t::ok;
	}

	//double frequency_;
	int16_t adding_coeff;
	int16_t sample_prep;
	std::mutex frequency_mutex_;
};

class TriangleWaveGenerator : public WaveGenerator
{
public:
	TriangleWaveGenerator(double frequency) :WaveGenerator(pAudioFilter(), frequency),
		sample_prep(0), direction(true)
	{
	}
	/*
	void set_frequency(double frequency) {
		std::unique_lock<std::mutex> lock(frequency_mutex_);

		frequency_ = frequency;
	}
	*/
private:
	error_type_t do_process(audio_buffer_t& buffer)
	{
		std::unique_lock<std::mutex> lock(frequency_mutex_);
		const audio_params_t& params = buffer.params;
		adding_coeff = 65536 / (convert_rate_to_int(params.rate) / frequency_);
		for (auto& sample : buffer.data) {
			sample = sample_prep;
			int16_t temp;
			if (direction) {
				temp = sample_prep + adding_coeff;
				if (sample_prep > temp) {
					direction = false;
					sample_prep -= adding_coeff;
				}
				else {
					sample_prep += adding_coeff;
				}
			}
			else {
				temp = sample_prep - adding_coeff;
				if (sample_prep < temp) {
					direction = true;
					sample_prep += adding_coeff;
				}
				else {
					sample_prep -= adding_coeff;
				}
			}
		}
		buffer.valid_samples = buffer.data.size();
		return error_type_t::ok;
	}

	//double frequency_;
	int16_t adding_coeff;
	int16_t sample_prep;
	boolean direction;
	std::mutex frequency_mutex_;
};

// --------------------------

// ------ EFFECTS ------ 

class VolumeChanger : public AudioFilter
{
public:
	VolumeChanger(const pAudioFilter& child, double maxVolume) :AudioFilter(child),
		maxVolume_(maxVolume),time_(0.0)
	{}
	void set_volume(double volume)
	{
		std::unique_lock<std::mutex> lock(maxVolume_mutex_);
		maxVolume_ = volume;
	}

private:
	error_type_t do_process(audio_buffer_t& buffer)
	{
		const double step = 1.0 / convert_rate_to_int(buffer.params.rate);
		
		for (auto& sample : buffer.data) {
			sample *= maxVolume_;
			time_ = time_ + step;
		}
		buffer.valid_samples = buffer.data.size();
		return error_type_t::ok;
	}

	double time_;
	double maxVolume_;
	std::mutex maxVolume_mutex_;
};

/**
* @brief Creates tremolo effect - changes volume of buffer based on sine shape
* @param frequency - frequency of a modulating wave of type double in bounds from 0.0 to 20.0
* @param dept - how much volume should be taken in the peak of the wave: 1.0 = sound goes down to 0 volume, 0.0 = volume doesn't change at all; type double in bounds 0.0-1.0
*/
class TremoloEffect : public AudioFilter {
public:
	TremoloEffect(const pAudioFilter& child, double frequency, double dept, bool state) : AudioFilter(child), frequency(frequency), time(0.0), state(state) {
		if (dept > 1.0) {
			dept = 1.0;
		}
		else {
			if (dept < 0.0) {
				dept = 0.0;
			}
		}
		this->dept = dept / 2.0;
	}

	void set_state(bool state) {
		this->state = state;
	}
	void set_frequency(bool frequency) {
		this->frequency = frequency;
	}
	void set_dept(bool dept) {
		this->dept = dept;
	}
private:
	bool state;
	double frequency;
	double time;
	double dept;
	error_type_t do_process(audio_buffer_t& buffer) {
		const audio_params_t& params = buffer.params;
		const double step = 1.0 / convert_rate_to_int(buffer.params.rate);
		//const double step = (1 / seconds) / convert_rate_to_int(params.rate);
		if (state) {
			for (auto& sample : buffer.data) {
				double temp = this->dept * std::sin(time * frequency * pi2 + pi2 / 4) + 1 - this->dept;
				/*
				if (temp < 0) {
					temp *= -1;
				}
				*/
				sample *= temp;
				time = time + step;
			}
		}
		return error_type_t::ok;
	}
};

/**
* 
*/
class CircularBuffer {
private:
	size_t start;
	size_t end;
	size_t size;
	bool start_condition;
	audio_sample_t buffer[250];
public:
	CircularBuffer(size_t size) : size(size), start(0), end(0) , start_condition(true) {}
	
	void push_back(audio_sample_t sample) {
		if (start_condition) {
			if (start == end) {
				buffer[start] = sample;
				end++;
			}
			else {
				buffer[end] = sample;
				start_condition = false;
			}
		}
		else {
			// Teoreticky to nemusím modulovat, jelikož ve chvíli, co tu hodnotu přičítám, tak už ji moduluji.
			if ((end + 1) % size == (start % size)) {
				// Mám naplněný buffer
				buffer[start] = sample;
				end = (end + 1) % size;
				start = (start + 1) % size;
			}
			else {
				if (start < end) {
					end++;
					buffer[end] = sample;
				}
			}
		}
	}

	audio_sample_t read_sample(size_t index) {
		// čte se od startu
		return buffer[(start + index) % size];
	}
};

/**
* @brief Creates chorus effect on incoming audio buffer using circular buffer
* @param rate Rate of modulating wave, accepts input between 0-20 of type double
* @param delay_level Strength of modulating wave, accepts input between 1-10 of type double
* @param dry_wet_ratio Mix between two effects expressed in percentege from 0-1 of type double
*/
class ChorusEffect : public AudioFilter {
public:
	ChorusEffect(const pAudioFilter& child, double rate, double delay_level, double dry_wet_ratio, bool state, size_t circular_buffer_size = 250) : AudioFilter(child),
		rate(rate), delay_level(delay_level), dry_wet_ratio(dry_wet_ratio), time(0.0), start(true), state(state),
		circular_buffer(CircularBuffer(circular_buffer_size)), circular_buffer_size(circular_buffer_size) {
		//TODO: optimalizovat velikost circular bufferu, takhle to spíše jenom odhaduji
		if (rate < 0.0) {
			this->rate = 0.0;
		}
		else {
			if (rate > 20.0) {
				this->rate = 20.0;
			}
		}
		if (delay_level < 0.0) {
			this->delay_level = 0.0;
		}
		else {
			if (delay_level >10.0) {
				this->delay_level = 10.0;
			}
		}
		if (dry_wet_ratio < 0.0) {
			this->dry_wet_ratio = 0.0;
		}
		else {
			if (dry_wet_ratio > 1.0) {
				this->dry_wet_ratio = 1.0;
			}
		}
	}
	void set_state(bool state) {
		this->state = state;
	}
	void set_rate(double rate) {
		this->rate = rate;
	}
	void set_delay_level(double delay_level) {
		this->delay_level = delay_level;
	}
	void set_dry_wet_ration(double dry_wet_ratio) {
		this->dry_wet_ratio = dry_wet_ratio;
	}
private:
	void add_chorus(audio_sample_t& sample, const double& step) {
		// Spočtu jak moc velký delay mám mít
		double sin = std::sin(time * rate * pi2);
		size_t delay = ((std::sin(time * rate * pi2) + 1) / 2) * delay_level*20 +20;
		delay = circular_buffer_size - delay;
		sample.left = sample.left*0.4 + (circular_buffer.read_sample(delay).left * 0.4 * dry_wet_ratio);
		sample.right = sample.right * 0.4 + (circular_buffer.read_sample(delay-25).right * 0.4 * dry_wet_ratio);
	}

	error_type_t do_process(audio_buffer_t& buffer) {
		// V první fázi čekám na naplnění circular bufferu
		// Až bude naplněný, tak spustím funkci 
		const audio_params_t& params = buffer.params;
		const double step = 1.0 / convert_rate_to_int(buffer.params.rate);
		size_t counter = 0;
		std::unique_lock<std::mutex> lock(mutex_);
		if (state) {
			for (auto& sample : buffer.data) {
				if (start) {
					circular_buffer.push_back(sample);
					if (counter < circular_buffer_size) {
						counter++;
					}
					else {
						start = false;
					}
				}
				else {
					circular_buffer.push_back(sample);
					add_chorus(sample, step);
					//sample = circular_buffer.read_sample(0);
				}
				time = time + step;
			}
		}
		buffer.valid_samples = buffer.data.size();
		return error_type_t::ok;
	}
	bool state;
	double rate;
	double delay_level;
	double dry_wet_ratio;
	size_t circular_buffer_size;
	CircularBuffer circular_buffer;
	std::mutex mutex_;
	double time;
	bool start;
}; 

/**
* @brief Creates chorus effect on incoming audio buffer using circular buffer
* @param rate Rate of modulating wave, accepts input between 0-10 of type double
* @param delay_level Strength of modulating wave, accepts input between 1-10 of type double
* @param dry_wet_ratio Mix between two effects expressed in percentege from 0-100
*/
class ChorusEffectOld : public AudioFilter {
public:
	ChorusEffectOld(const pAudioFilter& child, double rate, double delay_level, double dry_wet_ratio, size_t circular_buffer_size = 100) : AudioFilter(child),
		rate(rate), delay_level(delay_level), dry_wet_ratio(dry_wet_ratio), time(0.0),
		circular_buffer(circular_buffer_t<audio_sample_t>(circular_buffer_size)), circular_buffer_size(circular_buffer_size) {
		//TODO: optimalizovat velikost circular bufferu, takhle to spíše jenom odhaduji
	}
private:
	audio_sample_t get_sample_from_circ_buff(circular_buffer_t<audio_sample_t>& buffer, size_t& index) {
		if (buffer.start_ < buffer.end_) {
			return buffer.data[buffer.start_ + index];
		}
		else {
			// Vždycky tady bude buffer.start_ větší než buffer.end_
			return buffer.data[(buffer.start_ + index) % buffer.data.size()];
		}
	}

	void addChorus(audio_sample_t& sample, const double& step) {
		// Spočtu jak moc velký delay mám mít
		double sin = std::sin(time * rate * pi2);
		size_t delay = ((std::sin(time * rate * pi2) + 1) / 2) * delay_level + 15;
		delay = circular_buffer_size - delay;
		// Sample, který chci přičítat k tomu originalnimu
		audio_sample_t desired_sample = circular_buffer.data[((circular_buffer.start_ + delay) % circular_buffer.data.size())];
		sample = circular_buffer.start_;
	}

	error_type_t do_process(audio_buffer_t& buffer) {
		// V první fázi čekám na naplnění circular bufferu
		// Až bude naplněný, tak spustím funkci 

		const audio_params_t& params = buffer.params;
		const double step = 1.0 / convert_rate_to_int(buffer.params.rate);
		size_t counter = 0;
		audio_sample_t sample_from_circ_buff = audio_sample_t();
		for (auto& sample : buffer.data) {
			/*
			if (counter < 100) {
				circular_buffer.store_data(const_cast<const audio_sample_t*>(&sample), 1);
				counter++;
			}
			else {
				circular_buffer.get_data_block(&sample_from_circ_buff, 1); // tímhle vytahuji data z circular bufferu o velikosti jednoho bloku
				circular_buffer.store_data(const_cast<const audio_sample_t*>(&sample), 1);
				//sample = circular_buffer.data[(circular_buffer.start_+10) % circular_buffer_size];
				sample = sample_from_circ_buff;
				//addChorus(sample, step);

			}
			time = time + step;
			*/
		}
		//buffer.valid_samples = buffer.data.size();
		return error_type_t::ok;
	}
	double rate;
	double delay_level;
	double dry_wet_ratio;
	size_t circular_buffer_size;
	circular_buffer_t<audio_sample_t> circular_buffer;

	double time;
};

class FadeInM : public AudioFilter {
public:
	FadeInM(const pAudioFilter& child, double seconds) : AudioFilter(child), seconds(seconds), multiply_coeff(0) { 
	}
private:
	double seconds;
	double multiply_coeff;
	error_type_t do_process(audio_buffer_t& buffer) {
		const audio_params_t& params = buffer.params;
		const double step = (1 / seconds) / convert_rate_to_int(params.rate);
		for (auto& sample : buffer.data) {
			sample *= multiply_coeff;
			if (multiply_coeff < 1.0) {
				multiply_coeff += step;
			} else {
				if (multiply_coeff > 1.0) {
					multiply_coeff = 1.0;
				}
			}
		}
		//buffer.valid_samples = buffer.data.size();
		return error_type_t::ok;
	}
};

class FadeOutM : public AudioFilter {
public:
	FadeOutM(const pAudioFilter& child, double seconds) : AudioFilter(child), seconds(seconds), multiply_coeff(1) {
	}
private:
	double seconds;
	double multiply_coeff;
	error_type_t do_process(audio_buffer_t& buffer) {
		const audio_params_t& params = buffer.params;
		const double step = (1 / seconds) / convert_rate_to_int(params.rate);
		for (auto& sample : buffer.data) {
			sample *= multiply_coeff;
			if (multiply_coeff > 0.0) {
				multiply_coeff -= step;
			}
			else {
				if (multiply_coeff < 0.0) {
					multiply_coeff = 0.0;
				}
			}
		}
		//buffer.valid_samples = buffer.data.size();
		return error_type_t::ok;
	}
};

class FadeOut : public AudioFilter
{
public:
	FadeOut(const pAudioFilter& child) :AudioFilter(child)
	{};
private:
	error_type_t do_process(audio_buffer_t& buffer)
	{
		const double step = 1.0 / convert_rate_to_int(buffer.params.rate);
		auto sample_index_beg = int(buffer.data.size() * fade_out_beginning);
		auto fade_step = 1.0 / (buffer.data.size() - sample_index_beg);
		double volume_level = 1.0;

		for (size_t i = sample_index_beg; i < buffer.data.size(); ++i) {
			volume_level -= fade_step;
			buffer.data[i] *= volume_level;
		}
		buffer.valid_samples = buffer.data.size();

		return error_type_t::ok;
	}
double time_;
double fade_out_beginning = 0.9 ;
};

// ---------------------------------------

class Envelope_effect : public AudioFilter
{
public:
	Envelope_effect(const pAudioFilter& child) :AudioFilter(child)
	{};

	void set_input(std::vector<std::pair<double, double>>& line) {
		input = line;
	}

	// When envelope effect is turned off that means that we want to reset sound volume to 100%
	void turn_off() {

	}

private:
	error_type_t do_process(audio_buffer_t& buffer)
	{
		return error_type_t::ok;
	}

	double volume_rate_;
	std::mutex volume_rate_mutex_;
	std::vector<std::pair<double, double>> input;
};


class Control : public iimavlib::SDLDevice, public AudioFilter
{
public:
	static const rgb_t black;
	static const rgb_t light_green;
	static const rgb_t dark_green;
	static const rgb_t red;
	static const rgb_t light_violet;
	static const rgb_t dark_violet;
	static const rgb_t light_blue;
	static const rgb_t dark_blue;
	static const rgb_t light_yellow;
	static const rgb_t dark_yellow;

	Control(const pAudioFilter& child, int width, int height) :
		SDLDevice(width, height, "Frequency control", false),
		AudioFilter(child),
		data_(rectangle_t(0, 0, width, height), black),
		position_(0, 0)
	{
		rgb_t color(0, 0, 255);

		// Draw the buttons
		iimavlib::draw_empty_rectangle(data_, envelope_input_window, border_thickness, dark_green);
		iimavlib::draw_circle(data_, envelope_off_button, red);

		iimavlib::draw_empty_rectangle(data_, chorus_input_1, border_thickness, dark_violet);
		iimavlib::draw_empty_rectangle(data_, chorus_input_2, border_thickness, dark_violet);
		iimavlib::draw_empty_rectangle(data_, chorus_input_3, border_thickness, dark_violet);
		iimavlib::draw_circle(data_, chorus_off_button, red);

		iimavlib::draw_empty_rectangle(data_, tremolo_input_1, border_thickness, dark_blue);
		iimavlib::draw_empty_rectangle(data_, tremolo_input_2, border_thickness, dark_blue);
		iimavlib::draw_circle(data_, tremolo_off_button, red);

		blit(data_);

		// Start the rendering thread
		start();
	}
	~Control() {
		// Stop the rendering thread
		stop();
	}
private:
	/// Video data
	iimavlib::video_buffer_t data_;
	/// position of the mouse
	std::pair<int, int> position_;
	/// Mutex to lock @em index_ and @em position_
	std::mutex position_mutex_;
	bool envelope_input_window_entred = false;
	bool envelope_effect_enabled = true;
	std::vector<std::pair<int, int>> envelope_input_line = std::vector<std::pair<int, int>>();

	// UI
	const rectangle_t envelope_input_window = rectangle_t(37, 60, 300, 190);
	const int border_thickness = 2;
	const rectangle_t envelope_input_window_without_borders = rectangle_t(
		envelope_input_window.x + border_thickness,
		envelope_input_window.y + border_thickness,
		envelope_input_window.width - border_thickness * 2,
		envelope_input_window.height - border_thickness * 2);

	const rectangle_t chorus_input_1 = rectangle_t(430, 405, 285, 35);
	const rectangle_t chorus_input_2 = rectangle_t(430, 450, 285, 35);
	const rectangle_t chorus_input_3 = rectangle_t(430, 495, 285, 35);

	const rectangle_t tremolo_input_1 = rectangle_t(430, 110, 285, 35);
	const rectangle_t tremolo_input_2 = rectangle_t(430, 155, 285, 35);

	const rectangle_t envelope_off_button = rectangle_t(355, 213, 30, 30);
	const rectangle_t chorus_off_button = rectangle_t(738, 504, 30, 30);
	const rectangle_t tremolo_off_button = rectangle_t(738, 213, 30, 30);

	// Filter controllers
	void set_frequency(double freq) {
		auto generator = std::dynamic_pointer_cast<WaveShapesGenerator>(get_child(4));
		generator->set_frequency(freq);
	}

	void set_envelope_effect_input(std::vector<std::pair<double, double>>& line) {
		auto envelope_effect = std::dynamic_pointer_cast<Envelope_effect>(get_child(2));
		envelope_effect->set_input(line);
	}

	void turn_off_envelope_effect() {
		auto envelope_effect = std::dynamic_pointer_cast<Envelope_effect>(get_child(2));
		envelope_effect->turn_off();
	}

	void turn_off_tremolo() {

	}
	/**
	 * Overloaded method for handling keys from SDL window
	 * @param key  Number of the key pressed, defined in keys.h
	 * @param pressed True if the key was pressed, false if the key was released
	 * @return false if the program should end, true otherwise
	 */
	bool do_key_pressed(const int key, bool pressed) {
		using namespace iimavlib::keys;
		if (pressed) {
			auto note = notes.find(key);
			if (note != notes.end()) {
				set_frequency(note->second);
			}
			else {
				switch (key) {
					// If key Q or ESCAPE was pressed, we want to exit the program
				case 'q':
				case key_escape: return false;
				}
			}
		}
		update_screen(); // Update the screen　immediately to reflect the keypress
		return true;
	}

#define DBOUT( s )            \
{                             \
   std::ostringstream os_;    \
   os_ << s;                   \
   OutputDebugString( os_.str().c_str() );  \
}
	/**
	 * Overloaded method for processing mouse buttons.
	 * @param button Index of button that triggered the event
	 * @param pressed true if button was pressed, false if released
	 * @return true, unless some fatal error has occurred
	 */
	virtual bool do_mouse_button(const int button, const bool pressed, const int x, const int  y)
	{
		if (button == 0 && pressed) {
			{
				// Check if it is inside envelope effect input window
				if (has_intersect(x, y, envelope_input_window) && envelope_effect_enabled) {
					hide_last_drawn_line();

					envelope_input_line.push_back(std::make_pair(x, y));
					envelope_input_window_entred = true;
				}

				// Check if it is the envelope effect off button
				if (has_intersect(x, y, envelope_off_button)) {
					if (envelope_effect_enabled) {
						envelope_effect_enabled = false;
						hide_last_drawn_line();
						iimavlib::draw_circle(data_, envelope_off_button, red);
						iimavlib::draw_empty_rectangle(data_, envelope_input_window, border_thickness, dark_green);
					}
					else {
						envelope_effect_enabled = true;
						iimavlib::draw_circle(data_, envelope_off_button, light_green);
						iimavlib::draw_empty_rectangle(data_, envelope_input_window, border_thickness, light_green);
					}
				}

				std::unique_lock<std::mutex> lock(position_mutex_); // Lock the variables

			}

			update_screen();
		}

		if (button == 0 && !pressed) {
			envelope_input_window_entred = false;
			convert_line_to_envelope_effect_input();
			envelope_input_line.clear();
		}

		return true;
	}

	virtual bool do_mouse_moved(const int x, const int y, const int dx, const int dy) {
		if (envelope_input_window_entred && has_intersect(x, y, envelope_input_window)) {
			iimavlib::draw_line_thick(data_, rectangle_t(envelope_input_line.back().first,
				envelope_input_line.back().second, 0, 0),
				rectangle_t(x, y, 0, 0), 2, light_green);

			envelope_input_line.push_back(std::make_pair(x, y));
		}

		update_screen();
		return true;
	}

	void update_screen()
	{
		// And push it to the rendering thread
		blit(data_);
	}

	error_type_t do_process(iimavlib::audio_buffer_t& /*buffer*/)
	{
		if (is_stopped()) return error_type_t::failed;
		// Not touching the data, simply passing it through
		return error_type_t::ok;
	}

	bool has_intersect(int mouse_x, int mouse_y, rectangle_t active_zone) {
		if (mouse_x > active_zone.x && mouse_x < active_zone.width + active_zone.x &&
			mouse_y > active_zone.y && mouse_y < active_zone.height + active_zone.y) {
			return true;
		}
		return false;
	}

	void hide_last_drawn_line() {
		iimavlib::draw_rectangle(data_, envelope_input_window_without_borders, black);
	}

	void convert_line_to_envelope_effect_input() {
		// Pixel on x axis is 0.02 sec, from 0 sec to 5 sec on the whole axis
		double axis_x_rate = 5;
		// Pixel on y axis is 0.66 % of the max volume, from 0% to 100% on the whole axis

		std::vector<std::pair<double, double>> envelope_converted_input = std::vector<std::pair<double, double>>();
		for (auto p : envelope_input_line) {
			// Offset
			p.first -= envelope_input_window.x;
			p.second -= envelope_input_window.y;

			envelope_converted_input.push_back(std::make_pair(axis_x_rate * ((double)p.first / envelope_input_window.width),
				1 - ((double)p.second / envelope_input_window.height)));

		}

		set_envelope_effect_input(envelope_converted_input);
	}
};

const rgb_t Control::black(0, 0, 0);
const rgb_t Control::light_green(50, 205, 50);
const rgb_t Control::dark_green(0, 100, 0);
const rgb_t Control::red(255, 0, 0);
const rgb_t Control::light_violet(177, 25, 175);
const rgb_t Control::dark_violet(100, 30, 100);
const rgb_t Control::light_yellow(220, 200, 15);
const rgb_t Control::dark_yellow(115, 105, 30);
const rgb_t Control::light_blue(0, 165, 235);
const rgb_t Control::dark_blue(40, 105, 132);


int main(int argc, char** argv) try
{
	/* ******************************************************************
	 *                      Process parameters
	 ****************************************************************** */
	//if (argc < 2) {
	//	logger[log_level::fatal] << "Not enough parameters. Specify the frequency, please.";
	//	logger[log_level::fatal] << "Usage: " << argv[0] << " frequency [audio_device]";
	//	return 1;
	//}
	//const double frequency = std::stod(argv[1]);
	//const double frequency = simple_cast<double>(argv[1]);
	const double frequency = 660;
	const double maxVolumeLevel = 1.0;
	logger[log_level::debug] << "Generating sine with frequency " << frequency << "Hz.";

	audio_id_t device_id = PlatformDevice::default_device();
	if (argc > 2) {
		device_id = simple_cast<audio_id_t>(argv[2]);
	}
	logger[log_level::debug] << "Using audio device " << device_id;

	/* ******************************************************************
	 *                      Create and run the filter chain
	 ****************************************************************** */

	 // Create filter chain


	auto chain = filter_chain<WaveShapesGenerator>(frequency, WaveShape::Square)
		.add<WaveSink>("xx.wav")
		//.add<VolumeChanger>(maxVolumeLevel)
		.add<Envelope_effect>()
		//.add<AddNoise>()
		//.add<FadeOutM>(1.0)
		.add<TremoloEffect>(1.1, 0.8, false)
		//.add<SineMultiply>(7)
		//.add<SimpleEchoFilter>(0.01, 0.9)
		.add<ChorusEffect>(5, 3, 0.8, true)
		//.add<SineMultiply>(329.63)
		//.add<SineMultiply>(392.00)
		//.add<VolumeChanger>()
		.add<Control>(800, 600)
		.add<PlatformSink>(device_id) // moje reproduktory
		.sink();

	// Start the filters
	chain->run();

	/*
	 * Alternative syntax would be:
	 * auto sine = std::make_shared<SineGenerator>(frequency);
	 * auto sink = std::make_shared<PlatformSink>(sine, device_id);
	 * sink->run();
	 *
	 */


}
catch (std::exception& e)
{
	logger[log_level::fatal] << "ERROR: An error occurred during program run: " << e.what();
}
