#include <inttypes.h>
#include <math.h>

typedef __uint128_t pcg128_t;

#define PCG_128BIT_CONSTANT(high,low) ((((pcg128_t)high) << 64) + low)
#define PCG_DEFAULT_MULTIPLIER_128 PCG_128BIT_CONSTANT(2549297995355413924ULL,4865540595714422341ULL)
#define PCG_DEFAULT_INCREMENT_128  PCG_128BIT_CONSTANT(6364136223846793005ULL,1442695040888963407ULL)

struct pcg64_state {
  pcg128_t state;
  pcg128_t inc;
};

static inline void pcg_setseq_128_step_r(struct pcg64_state* rng) {
  rng->state = rng->state * PCG_DEFAULT_MULTIPLIER_128 + rng->inc;
}

static inline uint64_t pcg_rotr_64(uint64_t value, unsigned int rot) {
  return (value >> rot) | (value << ((- rot) & 63));
}

static inline uint64_t pcg_output(pcg128_t state) {
  return pcg_rotr_64(((uint64_t)(state >> 64u)) ^ (uint64_t)state, (unsigned int)(state >> 122u));
}

static inline pcg128_t pcg_advance_lcg_128(pcg128_t state, pcg128_t delta, pcg128_t cur_mult, pcg128_t cur_plus) {
  pcg128_t acc_mult = 1u;
  pcg128_t acc_plus = 0u;
  while (delta > 0) {
    if (delta & 1) {
      acc_mult *= cur_mult;
      acc_plus = acc_plus * cur_mult + cur_plus;
    }
    cur_plus = (cur_mult + 1) * cur_plus;
    cur_mult *= cur_mult;
    delta /= 2;
  }
  return acc_mult * state + acc_plus;
}

void pcg64_srandom(struct pcg64_state* rng, uint64_t initstate_h, uint64_t initstate_l, uint64_t initseq_h, uint64_t initseq_l) {
//void pcg64_srandom(struct pcg64_state* rng, pcg128_t initstate, pcg128_t initseq) {
  pcg128_t initstate = PCG_128BIT_CONSTANT(initstate_h, initstate_l);
  pcg128_t initseq   = PCG_128BIT_CONSTANT(initseq_h, initseq_l);

  rng->state = 0U;
  rng->inc = (initseq << 1u) | 1u;
  pcg_setseq_128_step_r(rng);
  rng->state += initstate;
  pcg_setseq_128_step_r(rng);
}

uint64_t pcg64_random(struct pcg64_state* rng) {
  pcg_setseq_128_step_r(rng);
  return pcg_output(rng->state);
}

void pcg64_random_n(struct pcg64_state* rng, int64_t n, uint64_t *x) {
	int i;
	for (i = 0; i < n; ++i) {
		x[i] = pcg64_random(rng);
	}
}

uint64_t pcg64_boundedrand(struct pcg64_state* rng, uint64_t bound) {
  uint64_t threshold = -bound % bound;
  for (;;) {
    uint64_t r = pcg64_random(rng);
    if (r >= threshold) return r % bound;
  }
}

void pcg64_boundedrand_n(struct pcg64_state* rng, uint64_t bound, int64_t n, uint64_t *x) {
	int i;
	for (i = 0; i < n; ++i) {
		x[i] = pcg64_boundedrand(rng, bound);
	}
}

/*
 * random_real: Generate a stream of bits uniformly at random and
 * interpret it as the fractional part of the binary expansion of a
 * number in [0, 1], 0.00001010011111010100...; then round it.
 */
double pcg64_random_real(struct pcg64_state* rng) {
	int exponent = -64;
	uint64_t significand;
	unsigned shift;

	/*
	 * Read zeros into the exponent until we hit a one; the rest
	 * will go into the significand.
	 */
	while ((significand = pcg64_random(rng)) == 0) {
		exponent -= 64;
		/*
		 * If the exponent falls below -1074 = emin + 1 - p,
		 * the exponent of the smallest subnormal, we are
		 * guaranteed the result will be rounded to zero.  This
		 * case is so unlikely it will happen in realistic
		 * terms only if pcg64_random is broken.
		 */
		if (exponent < -1074) return 0;
	}

	/*
	 * There is a 1 somewhere in significand, not necessarily in
	 * the most significant position.  If there are leading zeros,
	 * shift them into the exponent and refill the less-significant
	 * bits of the significand.  Can't predict one way or another
	 * whether there are leading zeros: there's a fifty-fifty
	 * chance, if pcg64_random is uniformly distributed.
	 */
  //_BitScanReverse64(&shift, significand);
  uint64_t bsr;
  asm ("bsrq %[significand], %[bsr]" : [bsr] "=r" (bsr) : [significand] "r" (significand) : "cc");
  shift = 63 - bsr;
	if (shift != 0) {
		exponent -= shift;
		significand <<= shift;
		significand |= (pcg64_random(rng) >> (64 - shift));
	}

	/*
	 * Set the sticky bit, since there is almost surely another 1
	 * in the bit stream.  Otherwise, we might round what looks
	 * like a tie to even when, almost surely, were we to look
	 * further in the bit stream, there would be a 1 breaking the
	 * tie.
	 */
	significand |= 1;

	/*
	 * Finally, convert to double (rounding) and scale by
	 * 2^exponent.
	 */
	return ldexp((double)significand, exponent);
}

void pcg64_random_real_n(struct pcg64_state* rng, int64_t n, double *x) {
	int i;
	for (i = 0; i < n; ++i) {
		x[i] = pcg64_random_real(rng);
	}
}

void pcg64_advance(struct pcg64_state* rng, uint64_t delta_h, uint64_t delta_l) {
//void pcg64_advance(struct pcg64_state* rng, pcg128_t delta) {
  pcg128_t delta = PCG_128BIT_CONSTANT(delta_h, delta_l);
  rng->state = pcg_advance_lcg_128(rng->state, delta, PCG_DEFAULT_MULTIPLIER_128, rng->inc);
}
