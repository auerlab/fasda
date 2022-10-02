// FIXME: Just a skeleton
void    fc_count(void);


/*
 *  Generate all combinations n choose 2 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge2(double fc_list[], unsigned long fc_count,
			double observed_fc_mean)

{
    unsigned long   fc_ge = 0;
    double          fc_mean;
    
    unsigned long  c1, c2;

    for (c1 = 0; c1 < fc_count; ++c1)
     for (c2 = c1 + 1; c2 < fc_count; ++c2)
     {
         fc_mean = (fc_list[c1] + fc_list[c2]) / 2;
         if ( fc_mean >= observed_fc_mean ) ++fc_ge;
     }
    return fc_ge;
}

/*
 *  Generate all combinations n choose 3 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge3(double fc_list[], unsigned long fc_count,
			double observed_fc_mean)

{
    unsigned long   fc_ge = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3;

    for (c1 = 0; c1 < fc_count; ++c1)
     for (c2 = c1 + 1; c2 < fc_count; ++c2)
      for (c3 = c2 + 1; c3 < fc_count; ++c3)
      {
          fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3]) / 3;
          if ( fc_mean >= observed_fc_mean ) ++fc_ge;
      }
    return fc_ge;
}

/*
 *  Generate all combinations n choose 4 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge4(double fc_list[], unsigned long fc_count,
			double observed_fc_mean)

{
    unsigned long   fc_ge = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3, c4;

    for (c1 = 0; c1 < fc_count; ++c1)
     for (c2 = c1 + 1; c2 < fc_count; ++c2)
      for (c3 = c2 + 1; c3 < fc_count; ++c3)
       for (c4 = c3 + 1; c4 < fc_count; ++c4)
       {
           fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3] + fc_list[c4]) / 4;
           if ( fc_mean >= observed_fc_mean ) ++fc_ge;
       }
    return fc_ge;
}

/*
 *  Generate all combinations n choose 5 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge5(double fc_list[], unsigned long fc_count,
			double observed_fc_mean)

{
    unsigned long   fc_ge = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3, c4, c5;

    for (c1 = 0; c1 < fc_count; ++c1)
     for (c2 = c1 + 1; c2 < fc_count; ++c2)
      for (c3 = c2 + 1; c3 < fc_count; ++c3)
       for (c4 = c3 + 1; c4 < fc_count; ++c4)
        for (c5 = c4 + 1; c5 < fc_count; ++c5)
        {
            fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3] + fc_list[c4] + fc_list[c5]) / 5;
            if ( fc_mean >= observed_fc_mean ) ++fc_ge;
        }
    return fc_ge;
}

/*
 *  Generate all combinations n choose 6 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge6(double fc_list[], unsigned long fc_count,
			double observed_fc_mean)

{
    unsigned long   fc_ge = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3, c4, c5, c6;

    for (c1 = 0; c1 < fc_count; ++c1)
     for (c2 = c1 + 1; c2 < fc_count; ++c2)
      for (c3 = c2 + 1; c3 < fc_count; ++c3)
       for (c4 = c3 + 1; c4 < fc_count; ++c4)
        for (c5 = c4 + 1; c5 < fc_count; ++c5)
         for (c6 = c5 + 1; c6 < fc_count; ++c6)
         {
             fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3] + fc_list[c4] + fc_list[c5] + fc_list[c6]) / 6;
             if ( fc_mean >= observed_fc_mean ) ++fc_ge;
         }
    return fc_ge;
}

/*
 *  Generate all combinations n choose 7 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge7(double fc_list[], unsigned long fc_count,
			double observed_fc_mean)

{
    unsigned long   fc_ge = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3, c4, c5, c6, c7;

    for (c1 = 0; c1 < fc_count; ++c1)
     for (c2 = c1 + 1; c2 < fc_count; ++c2)
      for (c3 = c2 + 1; c3 < fc_count; ++c3)
       for (c4 = c3 + 1; c4 < fc_count; ++c4)
        for (c5 = c4 + 1; c5 < fc_count; ++c5)
         for (c6 = c5 + 1; c6 < fc_count; ++c6)
          for (c7 = c6 + 1; c7 < fc_count; ++c7)
          {
              fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3] + fc_list[c4] + fc_list[c5] + fc_list[c6] + fc_list[c7]) / 7;
              if ( fc_mean >= observed_fc_mean ) ++fc_ge;
          }
    return fc_ge;
}

/*
 *  Generate all combinations n choose 8 and count FCs >= observed.
 *  This is much faster than generic algorithms for generating n choose
 *  k lists for any n and k.
 */

unsigned long   fc_ge8(double fc_list[], unsigned long fc_count,
			double observed_fc_mean)

{
    unsigned long   fc_ge = 0;
    double          fc_mean;
    
    unsigned long  c1, c2, c3, c4, c5, c6, c7, c8;

    for (c1 = 0; c1 < fc_count; ++c1)
     for (c2 = c1 + 1; c2 < fc_count; ++c2)
      for (c3 = c2 + 1; c3 < fc_count; ++c3)
       for (c4 = c3 + 1; c4 < fc_count; ++c4)
        for (c5 = c4 + 1; c5 < fc_count; ++c5)
         for (c6 = c5 + 1; c6 < fc_count; ++c6)
          for (c7 = c6 + 1; c7 < fc_count; ++c7)
           for (c8 = c7 + 1; c8 < fc_count; ++c8)
           {
               fc_mean = (fc_list[c1] + fc_list[c2] + fc_list[c3] + fc_list[c4] + fc_list[c5] + fc_list[c6] + fc_list[c7] + fc_list[c8]) / 8;
               if ( fc_mean >= observed_fc_mean ) ++fc_ge;
           }
    return fc_ge;
}


unsigned long   fc_ge_count(double fc_list[], unsigned long fc_count,
		      unsigned long replicates, double observed_fc_mean)

{
    static unsigned long (*fc_ge_funcs[])(double fc_list[],
				    unsigned long fc_count,
				    double observed_fc_mean) =
    {
        fc_ge2,
        fc_ge3,
        fc_ge4,
        fc_ge5,
        fc_ge6,
        fc_ge7,
        fc_ge8
    };
    unsigned long  func_index = replicates - 2;
    
    return fc_ge_funcs[func_index](fc_list, fc_count, observed_fc_mean);
}
