/*
 *  Generated by combgen.sh.  DO NOT EDIT.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>


/*
 *  Generate all combinations n choose 2 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge2(double fc_list[], unsigned long fc_count,
			double observed_fc_mean, unsigned long *fc_mean_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 1, count = 0;
    double          fc_mean;
    
    unsigned long  c1, c2;

    for (c1 = 0; c1 < fc_count; c1 += increment)
     for (c2 = c1 + 1; c2 < fc_count; c2 += increment)
     {
         fc_mean = (fc_list[c1] + fc_list[c2]) / 2;
         if ( fc_mean >= observed_fc_mean ) ++fc_ge;
         else if ( 1.0 / fc_mean >= observed_fc_mean ) ++fc_le;
         ++count;
     }
    printf("FC means > observed = %lu  FC means < 1/observed = %lu  Ratio = %f\n",
	    fc_ge, fc_le, (double)fc_ge / fc_le);
    *fc_mean_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 3 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge3(double fc_list[], unsigned long fc_count,
			double observed_fc_mean, unsigned long *fc_mean_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 1, count = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3;

    for (c1 = 0; c1 < fc_count; c1 += increment)
     for (c2 = c1 + 1; c2 < fc_count; c2 += increment)
      for (c3 = c2 + 1; c3 < fc_count; c3 += increment)
      {
          fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3]) / 3;
          if ( fc_mean >= observed_fc_mean ) ++fc_ge;
          else if ( 1.0 / fc_mean >= observed_fc_mean ) ++fc_le;
          ++count;
      }
    printf("FC means > observed = %lu  FC means < 1/observed = %lu  Ratio = %f\n",
	    fc_ge, fc_le, (double)fc_ge / fc_le);
    *fc_mean_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 4 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge4(double fc_list[], unsigned long fc_count,
			double observed_fc_mean, unsigned long *fc_mean_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 1, count = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3, c4;

    for (c1 = 0; c1 < fc_count; c1 += increment)
     for (c2 = c1 + 1; c2 < fc_count; c2 += increment)
      for (c3 = c2 + 1; c3 < fc_count; c3 += increment)
       for (c4 = c3 + 1; c4 < fc_count; c4 += increment)
       {
           fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3] + fc_list[c4]) / 4;
           if ( fc_mean >= observed_fc_mean ) ++fc_ge;
           else if ( 1.0 / fc_mean >= observed_fc_mean ) ++fc_le;
           ++count;
       }
    printf("FC means > observed = %lu  FC means < 1/observed = %lu  Ratio = %f\n",
	    fc_ge, fc_le, (double)fc_ge / fc_le);
    *fc_mean_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 5 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge5(double fc_list[], unsigned long fc_count,
			double observed_fc_mean, unsigned long *fc_mean_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 1, count = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3, c4, c5;

    for (c1 = 0; c1 < fc_count; c1 += increment)
     for (c2 = c1 + 1; c2 < fc_count; c2 += increment)
      for (c3 = c2 + 1; c3 < fc_count; c3 += increment)
       for (c4 = c3 + 1; c4 < fc_count; c4 += increment)
        for (c5 = c4 + 1; c5 < fc_count; c5 += increment)
        {
            fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3] + fc_list[c4] + fc_list[c5]) / 5;
            if ( fc_mean >= observed_fc_mean ) ++fc_ge;
            else if ( 1.0 / fc_mean >= observed_fc_mean ) ++fc_le;
            ++count;
        }
    printf("FC means > observed = %lu  FC means < 1/observed = %lu  Ratio = %f\n",
	    fc_ge, fc_le, (double)fc_ge / fc_le);
    *fc_mean_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 6 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge6(double fc_list[], unsigned long fc_count,
			double observed_fc_mean, unsigned long *fc_mean_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 2, count = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3, c4, c5, c6;

    for (c1 = 0; c1 < fc_count; c1 += increment)
     for (c2 = c1 + 1 + random() % increment; c2 < fc_count; c2 += increment)
      for (c3 = c2 + 1 + random() % increment; c3 < fc_count; c3 += increment)
       for (c4 = c3 + 1 + random() % increment; c4 < fc_count; c4 += increment)
        for (c5 = c4 + 1 + random() % increment; c5 < fc_count; c5 += increment)
         for (c6 = c5 + 1 + random() % increment; c6 < fc_count; c6 += increment)
         {
             fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3] + fc_list[c4] + fc_list[c5] + fc_list[c6]) / 6;
             if ( fc_mean >= observed_fc_mean ) ++fc_ge;
             else if ( 1.0 / fc_mean >= observed_fc_mean ) ++fc_le;
             ++count;
         }
    printf("FC means > observed = %lu  FC means < 1/observed = %lu  Ratio = %f\n",
	    fc_ge, fc_le, (double)fc_ge / fc_le);
    *fc_mean_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 7 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge7(double fc_list[], unsigned long fc_count,
			double observed_fc_mean, unsigned long *fc_mean_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 4, count = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3, c4, c5, c6, c7;

    for (c1 = 0; c1 < fc_count; c1 += increment)
     for (c2 = c1 + 1 + random() % increment; c2 < fc_count; c2 += increment)
      for (c3 = c2 + 1 + random() % increment; c3 < fc_count; c3 += increment)
       for (c4 = c3 + 1 + random() % increment; c4 < fc_count; c4 += increment)
        for (c5 = c4 + 1 + random() % increment; c5 < fc_count; c5 += increment)
         for (c6 = c5 + 1 + random() % increment; c6 < fc_count; c6 += increment)
          for (c7 = c6 + 1 + random() % increment; c7 < fc_count; c7 += increment)
          {
              fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3] + fc_list[c4] + fc_list[c5] + fc_list[c6] + fc_list[c7]) / 7;
              if ( fc_mean >= observed_fc_mean ) ++fc_ge;
              else if ( 1.0 / fc_mean >= observed_fc_mean ) ++fc_le;
              ++count;
          }
    printf("FC means > observed = %lu  FC means < 1/observed = %lu  Ratio = %f\n",
	    fc_ge, fc_le, (double)fc_ge / fc_le);
    *fc_mean_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 8 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge8(double fc_list[], unsigned long fc_count,
			double observed_fc_mean, unsigned long *fc_mean_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 6, count = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3, c4, c5, c6, c7, c8;

    for (c1 = 0; c1 < fc_count; c1 += increment)
     for (c2 = c1 + 1 + random() % increment; c2 < fc_count; c2 += increment)
      for (c3 = c2 + 1 + random() % increment; c3 < fc_count; c3 += increment)
       for (c4 = c3 + 1 + random() % increment; c4 < fc_count; c4 += increment)
        for (c5 = c4 + 1 + random() % increment; c5 < fc_count; c5 += increment)
         for (c6 = c5 + 1 + random() % increment; c6 < fc_count; c6 += increment)
          for (c7 = c6 + 1 + random() % increment; c7 < fc_count; c7 += increment)
           for (c8 = c7 + 1 + random() % increment; c8 < fc_count; c8 += increment)
           {
               fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3] + fc_list[c4] + fc_list[c5] + fc_list[c6] + fc_list[c7] + fc_list[c8]) / 8;
               if ( fc_mean >= observed_fc_mean ) ++fc_ge;
               else if ( 1.0 / fc_mean >= observed_fc_mean ) ++fc_le;
               ++count;
           }
    printf("FC means > observed = %lu  FC means < 1/observed = %lu  Ratio = %f\n",
	    fc_ge, fc_le, (double)fc_ge / fc_le);
    *fc_mean_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 9 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge9(double fc_list[], unsigned long fc_count,
			double observed_fc_mean, unsigned long *fc_mean_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 10, count = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3, c4, c5, c6, c7, c8, c9;

    for (c1 = 0; c1 < fc_count; c1 += increment)
     for (c2 = c1 + 1 + random() % increment; c2 < fc_count; c2 += increment)
      for (c3 = c2 + 1 + random() % increment; c3 < fc_count; c3 += increment)
       for (c4 = c3 + 1 + random() % increment; c4 < fc_count; c4 += increment)
        for (c5 = c4 + 1 + random() % increment; c5 < fc_count; c5 += increment)
         for (c6 = c5 + 1 + random() % increment; c6 < fc_count; c6 += increment)
          for (c7 = c6 + 1 + random() % increment; c7 < fc_count; c7 += increment)
           for (c8 = c7 + 1 + random() % increment; c8 < fc_count; c8 += increment)
            for (c9 = c8 + 1 + random() % increment; c9 < fc_count; c9 += increment)
            {
                fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3] + fc_list[c4] + fc_list[c5] + fc_list[c6] + fc_list[c7] + fc_list[c8] + fc_list[c9]) / 9;
                if ( fc_mean >= observed_fc_mean ) ++fc_ge;
                else if ( 1.0 / fc_mean >= observed_fc_mean ) ++fc_le;
                ++count;
            }
    printf("FC means > observed = %lu  FC means < 1/observed = %lu  Ratio = %f\n",
	    fc_ge, fc_le, (double)fc_ge / fc_le);
    *fc_mean_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 10 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge10(double fc_list[], unsigned long fc_count,
			double observed_fc_mean, unsigned long *fc_mean_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 15, count = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3, c4, c5, c6, c7, c8, c9, c10;

    for (c1 = 0; c1 < fc_count; c1 += increment)
     for (c2 = c1 + 1 + random() % increment; c2 < fc_count; c2 += increment)
      for (c3 = c2 + 1 + random() % increment; c3 < fc_count; c3 += increment)
       for (c4 = c3 + 1 + random() % increment; c4 < fc_count; c4 += increment)
        for (c5 = c4 + 1 + random() % increment; c5 < fc_count; c5 += increment)
         for (c6 = c5 + 1 + random() % increment; c6 < fc_count; c6 += increment)
          for (c7 = c6 + 1 + random() % increment; c7 < fc_count; c7 += increment)
           for (c8 = c7 + 1 + random() % increment; c8 < fc_count; c8 += increment)
            for (c9 = c8 + 1 + random() % increment; c9 < fc_count; c9 += increment)
             for (c10 = c9 + 1 + random() % increment; c10 < fc_count; c10 += increment)
             {
                 fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3] + fc_list[c4] + fc_list[c5] + fc_list[c6] + fc_list[c7] + fc_list[c8] + fc_list[c9] + fc_list[c10]) / 10;
                 if ( fc_mean >= observed_fc_mean ) ++fc_ge;
                 else if ( 1.0 / fc_mean >= observed_fc_mean ) ++fc_le;
                 ++count;
             }
    printf("FC means > observed = %lu  FC means < 1/observed = %lu  Ratio = %f\n",
	    fc_ge, fc_le, (double)fc_ge / fc_le);
    *fc_mean_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}


unsigned long   fc_ge_count(double fc_list[], unsigned long fc_count,
		      unsigned long replicates, double observed_fc_mean,
		      unsigned long *fc_mean_count)

{
    static unsigned long (*fc_ge_funcs[])(double fc_list[],
				    unsigned long fc_count,
				    double observed_fc_mean,
				    unsigned long *fc_mean_count) =
    {
        fc_ge2,
        fc_ge3,
        fc_ge4,
        fc_ge5,
        fc_ge6,
        fc_ge7,
        fc_ge8,
        fc_ge9,
        fc_ge10
    };
    unsigned long  func_index = replicates - 2;
    
    srandom(time(NULL));
    return fc_ge_funcs[func_index](fc_list, fc_count,
				   observed_fc_mean, fc_mean_count);
}