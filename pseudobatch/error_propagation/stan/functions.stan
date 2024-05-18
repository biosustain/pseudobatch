/*  Cumulative sum of a vector.

    Use this function to calculate accumulated dilution factor.
*/
vector cumulative_product_vector(vector v){
  vector[rows(v)] out;
  real cumulative_product = 1;
  for (i in 1:rows(v)){
    cumulative_product *= v[i];
    out[i] = cumulative_product;
  }
  return out;
}

/* Function that does the pseudo-batch transformation */
vector pseudobatch_transform(vector v,   // before sampling volumes
                             vector s,   // sample volumes
                             vector c,   // target concentration
                             vector f,   // feed in previous interval
                             real cfeed  // concentration of feed
                             ) {
  // After sample reactor volume is before sample volume minus sample volume.
  vector[rows(v)] a = v - s;
  // Dilution factor is s + a divided by a shifted once (first entry filled by 1).
  vector[rows(v)] df = append_row(1, (s[2:] + a[2:]) ./ a[:rows(v)-1]);
  // Accumulated dilution factor is cumulative product of dilution factors.
  vector[rows(v)] adf = cumulative_product_vector(df);
  // Output is target concentration times adf times cumulative sume of feed
  // concentration times adf times accumlated feed volume divided by total
  // (before sampling) volume.
  return c .* adf - cumulative_sum(cfeed * adf .* f ./ v);
}

int count_zeros(vector v){
  int out = 0;
  for (vv in v){
    if (vv == 0) out += 1;
  }
  return out;
}

array[] int vector_zero_ixs(vector v){
  /* Find the indexes of v that are zero. */
  int N_out = count_zeros(v);
  array[N_out] int out;
  int iout = 1;
  for (iv in 1:size(v)){
    if (v[iv] == 0){
      out[iout] = iv;
      iout += 1;
    }
  }
  return out;
}

array[] int vector_nonzero_ixs(vector v){
  /* Find the indexes of v that are not zero. */
  int N_out = size(v) - count_zeros(v);
  array[N_out] int out;
  int iout = 1;
  for (iv in 1:size(v)){
    if (v[iv] != 0){
      out[iout] = iv;
      iout += 1;
    }
  }
  return out;
}
