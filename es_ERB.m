function ERB(fc,fm)

bw = 24.7 * (.00437*fc+1);
bw = bw * 1.25;
upper = fc + bw;
lower = fc - bw;
mod_upper = fc + fm;
mod_lower = fc - fm;
upper_overlap = upper - mod_upper
lower_overlap = mod_lower - lower