/*
 *  Generated by combgen.sh.  DO NOT EDIT.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "pval.h"


/*
 *  Generate all combinations n choose 2 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge2(count_pair_t count_pairs[], unsigned long pair_count,
			double observed_fc,
			unsigned long *fc_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 1, count = 0,
		    fc_g1 = 0, fc_l1 = 0;
    double          fc;
    
    unsigned long  c1, c2;

    for (c1 = 0; c1 < pair_count; c1 += increment)
     for (c2 = c1 + 1; c2 < pair_count; c2 += increment)
     {
         fc = (double)(
                  count_pairs[c1].c2_count +
                  count_pairs[c2].c2_count
              )
              / 
              (
                  count_pairs[c1].c1_count +
                  count_pairs[c2].c1_count
              );
         if ( fc >= observed_fc ) ++fc_ge;
         else if ( fc <= 1.0 / observed_fc ) ++fc_le;
         if ( fc > 1 ) ++fc_g1;
         else if ( fc < 1 ) ++fc_l1;
         ++count;
     }
    printf("FCs > 1 = %-5lu           FCs < 1 = %-5lu\n", fc_g1, fc_l1);
    printf("FCs > %0.5f = %-5lu     FCs < %0.5f = %-5lu\n",
	    observed_fc, fc_ge, 1.0 / observed_fc, fc_le);
    *fc_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 3 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge3(count_pair_t count_pairs[], unsigned long pair_count,
			double observed_fc,
			unsigned long *fc_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 1, count = 0,
		    fc_g1 = 0, fc_l1 = 0;
    double          fc;
    
    unsigned long  c1, c2, c3;

    for (c1 = 0; c1 < pair_count; c1 += increment)
     for (c2 = c1 + 1; c2 < pair_count; c2 += increment)
      for (c3 = c2 + 1; c3 < pair_count; c3 += increment)
      {
          fc = (double)(
                   count_pairs[c1].c2_count +
                   count_pairs[c2].c2_count +
                   count_pairs[c3].c2_count
               )
               / 
               (
                   count_pairs[c1].c1_count +
                   count_pairs[c2].c1_count +
                  count_pairs[c3].c1_count
               );
          if ( fc >= observed_fc ) ++fc_ge;
          else if ( fc <= 1.0 / observed_fc ) ++fc_le;
          if ( fc > 1 ) ++fc_g1;
          else if ( fc < 1 ) ++fc_l1;
          ++count;
      }
    printf("FCs > 1 = %-5lu           FCs < 1 = %-5lu\n", fc_g1, fc_l1);
    printf("FCs > %0.5f = %-5lu     FCs < %0.5f = %-5lu\n",
	    observed_fc, fc_ge, 1.0 / observed_fc, fc_le);
    *fc_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 4 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge4(count_pair_t count_pairs[], unsigned long pair_count,
			double observed_fc,
			unsigned long *fc_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 1, count = 0,
		    fc_g1 = 0, fc_l1 = 0;
    double          fc;
    
    unsigned long  c1, c2, c3, c4;

    for (c1 = 0; c1 < pair_count; c1 += increment)
     for (c2 = c1 + 1; c2 < pair_count; c2 += increment)
      for (c3 = c2 + 1; c3 < pair_count; c3 += increment)
       for (c4 = c3 + 1; c4 < pair_count; c4 += increment)
       {
           fc = (double)(
                    count_pairs[c1].c2_count +
                    count_pairs[c2].c2_count +
                    count_pairs[c3].c2_count +
                    count_pairs[c4].c2_count
                )
                / 
                (
                    count_pairs[c1].c1_count +
                    count_pairs[c2].c1_count +
                    count_pairs[c3].c1_count +
                  count_pairs[c4].c1_count
                );
           if ( fc >= observed_fc ) ++fc_ge;
           else if ( fc <= 1.0 / observed_fc ) ++fc_le;
           if ( fc > 1 ) ++fc_g1;
           else if ( fc < 1 ) ++fc_l1;
           ++count;
       }
    printf("FCs > 1 = %-5lu           FCs < 1 = %-5lu\n", fc_g1, fc_l1);
    printf("FCs > %0.5f = %-5lu     FCs < %0.5f = %-5lu\n",
	    observed_fc, fc_ge, 1.0 / observed_fc, fc_le);
    *fc_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 5 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge5(count_pair_t count_pairs[], unsigned long pair_count,
			double observed_fc,
			unsigned long *fc_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 1, count = 0,
		    fc_g1 = 0, fc_l1 = 0;
    double          fc;
    
    unsigned long  c1, c2, c3, c4, c5;

    for (c1 = 0; c1 < pair_count; c1 += increment)
     for (c2 = c1 + 1; c2 < pair_count; c2 += increment)
      for (c3 = c2 + 1; c3 < pair_count; c3 += increment)
       for (c4 = c3 + 1; c4 < pair_count; c4 += increment)
        for (c5 = c4 + 1; c5 < pair_count; c5 += increment)
        {
            fc = (double)(
                     count_pairs[c1].c2_count +
                     count_pairs[c2].c2_count +
                     count_pairs[c3].c2_count +
                     count_pairs[c4].c2_count +
                     count_pairs[c5].c2_count
                 )
                 / 
                 (
                     count_pairs[c1].c1_count +
                     count_pairs[c2].c1_count +
                     count_pairs[c3].c1_count +
                     count_pairs[c4].c1_count +
                  count_pairs[c5].c1_count
                 );
            if ( fc >= observed_fc ) ++fc_ge;
            else if ( fc <= 1.0 / observed_fc ) ++fc_le;
            if ( fc > 1 ) ++fc_g1;
            else if ( fc < 1 ) ++fc_l1;
            ++count;
        }
    printf("FCs > 1 = %-5lu           FCs < 1 = %-5lu\n", fc_g1, fc_l1);
    printf("FCs > %0.5f = %-5lu     FCs < %0.5f = %-5lu\n",
	    observed_fc, fc_ge, 1.0 / observed_fc, fc_le);
    *fc_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 6 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge6(count_pair_t count_pairs[], unsigned long pair_count,
			double observed_fc,
			unsigned long *fc_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 2, count = 0,
		    fc_g1 = 0, fc_l1 = 0;
    double          fc;
    
    unsigned long  c1, c2, c3, c4, c5, c6;

    for (c1 = 0; c1 < pair_count; c1 += increment)
     for (c2 = c1 + 1 + random() % increment; c2 < pair_count; c2 += increment)
      for (c3 = c2 + 1 + random() % increment; c3 < pair_count; c3 += increment)
       for (c4 = c3 + 1 + random() % increment; c4 < pair_count; c4 += increment)
        for (c5 = c4 + 1 + random() % increment; c5 < pair_count; c5 += increment)
         for (c6 = c5 + 1 + random() % increment; c6 < pair_count; c6 += increment)
         {
             fc = (double)(
                      count_pairs[c1].c2_count +
                      count_pairs[c2].c2_count +
                      count_pairs[c3].c2_count +
                      count_pairs[c4].c2_count +
                      count_pairs[c5].c2_count +
                      count_pairs[c6].c2_count
                  )
                  / 
                  (
                      count_pairs[c1].c1_count +
                      count_pairs[c2].c1_count +
                      count_pairs[c3].c1_count +
                      count_pairs[c4].c1_count +
                      count_pairs[c5].c1_count +
                  count_pairs[c6].c1_count
                  );
             if ( fc >= observed_fc ) ++fc_ge;
             else if ( fc <= 1.0 / observed_fc ) ++fc_le;
             if ( fc > 1 ) ++fc_g1;
             else if ( fc < 1 ) ++fc_l1;
             ++count;
         }
    printf("FCs > 1 = %-5lu           FCs < 1 = %-5lu\n", fc_g1, fc_l1);
    printf("FCs > %0.5f = %-5lu     FCs < %0.5f = %-5lu\n",
	    observed_fc, fc_ge, 1.0 / observed_fc, fc_le);
    *fc_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 7 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge7(count_pair_t count_pairs[], unsigned long pair_count,
			double observed_fc,
			unsigned long *fc_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 4, count = 0,
		    fc_g1 = 0, fc_l1 = 0;
    double          fc;
    
    unsigned long  c1, c2, c3, c4, c5, c6, c7;

    for (c1 = 0; c1 < pair_count; c1 += increment)
     for (c2 = c1 + 1 + random() % increment; c2 < pair_count; c2 += increment)
      for (c3 = c2 + 1 + random() % increment; c3 < pair_count; c3 += increment)
       for (c4 = c3 + 1 + random() % increment; c4 < pair_count; c4 += increment)
        for (c5 = c4 + 1 + random() % increment; c5 < pair_count; c5 += increment)
         for (c6 = c5 + 1 + random() % increment; c6 < pair_count; c6 += increment)
          for (c7 = c6 + 1 + random() % increment; c7 < pair_count; c7 += increment)
          {
              fc = (double)(
                       count_pairs[c1].c2_count +
                       count_pairs[c2].c2_count +
                       count_pairs[c3].c2_count +
                       count_pairs[c4].c2_count +
                       count_pairs[c5].c2_count +
                       count_pairs[c6].c2_count +
                       count_pairs[c7].c2_count
                   )
                   / 
                   (
                       count_pairs[c1].c1_count +
                       count_pairs[c2].c1_count +
                       count_pairs[c3].c1_count +
                       count_pairs[c4].c1_count +
                       count_pairs[c5].c1_count +
                       count_pairs[c6].c1_count +
                  count_pairs[c7].c1_count
                   );
              if ( fc >= observed_fc ) ++fc_ge;
              else if ( fc <= 1.0 / observed_fc ) ++fc_le;
              if ( fc > 1 ) ++fc_g1;
              else if ( fc < 1 ) ++fc_l1;
              ++count;
          }
    printf("FCs > 1 = %-5lu           FCs < 1 = %-5lu\n", fc_g1, fc_l1);
    printf("FCs > %0.5f = %-5lu     FCs < %0.5f = %-5lu\n",
	    observed_fc, fc_ge, 1.0 / observed_fc, fc_le);
    *fc_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 8 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge8(count_pair_t count_pairs[], unsigned long pair_count,
			double observed_fc,
			unsigned long *fc_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 6, count = 0,
		    fc_g1 = 0, fc_l1 = 0;
    double          fc;
    
    unsigned long  c1, c2, c3, c4, c5, c6, c7, c8;

    for (c1 = 0; c1 < pair_count; c1 += increment)
     for (c2 = c1 + 1 + random() % increment; c2 < pair_count; c2 += increment)
      for (c3 = c2 + 1 + random() % increment; c3 < pair_count; c3 += increment)
       for (c4 = c3 + 1 + random() % increment; c4 < pair_count; c4 += increment)
        for (c5 = c4 + 1 + random() % increment; c5 < pair_count; c5 += increment)
         for (c6 = c5 + 1 + random() % increment; c6 < pair_count; c6 += increment)
          for (c7 = c6 + 1 + random() % increment; c7 < pair_count; c7 += increment)
           for (c8 = c7 + 1 + random() % increment; c8 < pair_count; c8 += increment)
           {
               fc = (double)(
                        count_pairs[c1].c2_count +
                        count_pairs[c2].c2_count +
                        count_pairs[c3].c2_count +
                        count_pairs[c4].c2_count +
                        count_pairs[c5].c2_count +
                        count_pairs[c6].c2_count +
                        count_pairs[c7].c2_count +
                        count_pairs[c8].c2_count
                    )
                    / 
                    (
                        count_pairs[c1].c1_count +
                        count_pairs[c2].c1_count +
                        count_pairs[c3].c1_count +
                        count_pairs[c4].c1_count +
                        count_pairs[c5].c1_count +
                        count_pairs[c6].c1_count +
                        count_pairs[c7].c1_count +
                  count_pairs[c8].c1_count
                    );
               if ( fc >= observed_fc ) ++fc_ge;
               else if ( fc <= 1.0 / observed_fc ) ++fc_le;
               if ( fc > 1 ) ++fc_g1;
               else if ( fc < 1 ) ++fc_l1;
               ++count;
           }
    printf("FCs > 1 = %-5lu           FCs < 1 = %-5lu\n", fc_g1, fc_l1);
    printf("FCs > %0.5f = %-5lu     FCs < %0.5f = %-5lu\n",
	    observed_fc, fc_ge, 1.0 / observed_fc, fc_le);
    *fc_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 9 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge9(count_pair_t count_pairs[], unsigned long pair_count,
			double observed_fc,
			unsigned long *fc_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 10, count = 0,
		    fc_g1 = 0, fc_l1 = 0;
    double          fc;
    
    unsigned long  c1, c2, c3, c4, c5, c6, c7, c8, c9;

    for (c1 = 0; c1 < pair_count; c1 += increment)
     for (c2 = c1 + 1 + random() % increment; c2 < pair_count; c2 += increment)
      for (c3 = c2 + 1 + random() % increment; c3 < pair_count; c3 += increment)
       for (c4 = c3 + 1 + random() % increment; c4 < pair_count; c4 += increment)
        for (c5 = c4 + 1 + random() % increment; c5 < pair_count; c5 += increment)
         for (c6 = c5 + 1 + random() % increment; c6 < pair_count; c6 += increment)
          for (c7 = c6 + 1 + random() % increment; c7 < pair_count; c7 += increment)
           for (c8 = c7 + 1 + random() % increment; c8 < pair_count; c8 += increment)
            for (c9 = c8 + 1 + random() % increment; c9 < pair_count; c9 += increment)
            {
                fc = (double)(
                         count_pairs[c1].c2_count +
                         count_pairs[c2].c2_count +
                         count_pairs[c3].c2_count +
                         count_pairs[c4].c2_count +
                         count_pairs[c5].c2_count +
                         count_pairs[c6].c2_count +
                         count_pairs[c7].c2_count +
                         count_pairs[c8].c2_count +
                         count_pairs[c9].c2_count
                     )
                     / 
                     (
                         count_pairs[c1].c1_count +
                         count_pairs[c2].c1_count +
                         count_pairs[c3].c1_count +
                         count_pairs[c4].c1_count +
                         count_pairs[c5].c1_count +
                         count_pairs[c6].c1_count +
                         count_pairs[c7].c1_count +
                         count_pairs[c8].c1_count +
                  count_pairs[c9].c1_count
                     );
                if ( fc >= observed_fc ) ++fc_ge;
                else if ( fc <= 1.0 / observed_fc ) ++fc_le;
                if ( fc > 1 ) ++fc_g1;
                else if ( fc < 1 ) ++fc_l1;
                ++count;
            }
    printf("FCs > 1 = %-5lu           FCs < 1 = %-5lu\n", fc_g1, fc_l1);
    printf("FCs > %0.5f = %-5lu     FCs < %0.5f = %-5lu\n",
	    observed_fc, fc_ge, 1.0 / observed_fc, fc_le);
    *fc_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}

/*
 *  Generate all combinations n choose 10 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge10(count_pair_t count_pairs[], unsigned long pair_count,
			double observed_fc,
			unsigned long *fc_count)

{
    // Using sample++ % sample_rate doesn't produce much gain
    // Go after loop increments instead
    unsigned long   fc_ge = 0, fc_le = 0, increment = 15, count = 0,
		    fc_g1 = 0, fc_l1 = 0;
    double          fc;
    
    unsigned long  c1, c2, c3, c4, c5, c6, c7, c8, c9, c10;

    for (c1 = 0; c1 < pair_count; c1 += increment)
     for (c2 = c1 + 1 + random() % increment; c2 < pair_count; c2 += increment)
      for (c3 = c2 + 1 + random() % increment; c3 < pair_count; c3 += increment)
       for (c4 = c3 + 1 + random() % increment; c4 < pair_count; c4 += increment)
        for (c5 = c4 + 1 + random() % increment; c5 < pair_count; c5 += increment)
         for (c6 = c5 + 1 + random() % increment; c6 < pair_count; c6 += increment)
          for (c7 = c6 + 1 + random() % increment; c7 < pair_count; c7 += increment)
           for (c8 = c7 + 1 + random() % increment; c8 < pair_count; c8 += increment)
            for (c9 = c8 + 1 + random() % increment; c9 < pair_count; c9 += increment)
             for (c10 = c9 + 1 + random() % increment; c10 < pair_count; c10 += increment)
             {
                 fc = (double)(
                          count_pairs[c1].c2_count +
                          count_pairs[c2].c2_count +
                          count_pairs[c3].c2_count +
                          count_pairs[c4].c2_count +
                          count_pairs[c5].c2_count +
                          count_pairs[c6].c2_count +
                          count_pairs[c7].c2_count +
                          count_pairs[c8].c2_count +
                          count_pairs[c9].c2_count +
                          count_pairs[c10].c2_count
                      )
                      / 
                      (
                          count_pairs[c1].c1_count +
                          count_pairs[c2].c1_count +
                          count_pairs[c3].c1_count +
                          count_pairs[c4].c1_count +
                          count_pairs[c5].c1_count +
                          count_pairs[c6].c1_count +
                          count_pairs[c7].c1_count +
                          count_pairs[c8].c1_count +
                          count_pairs[c9].c1_count +
                  count_pairs[c10].c1_count
                      );
                 if ( fc >= observed_fc ) ++fc_ge;
                 else if ( fc <= 1.0 / observed_fc ) ++fc_le;
                 if ( fc > 1 ) ++fc_g1;
                 else if ( fc < 1 ) ++fc_l1;
                 ++count;
             }
    printf("FCs > 1 = %-5lu           FCs < 1 = %-5lu\n", fc_g1, fc_l1);
    printf("FCs > %0.5f = %-5lu     FCs < %0.5f = %-5lu\n",
	    observed_fc, fc_ge, 1.0 / observed_fc, fc_le);
    *fc_count = count;
    
    // FIXME: Is this correct?
    return fc_ge + fc_le;
}


unsigned long   fc_ge_count(count_pair_t count_pairs[], unsigned long pair_count,
		      unsigned long replicates, double observed_fc,
		      unsigned long *fc_count)

{
    static unsigned long (*fc_ge_funcs[])(count_pair_t count_pairs[],
				    unsigned long pair_count,
				    double observed_fc,
				    unsigned long *fc_count) =
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
    return fc_ge_funcs[func_index](count_pairs, pair_count,
				   observed_fc, fc_count);
}
