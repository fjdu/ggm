from multiprocessing import Pool

def f(cmd):
    import os
    print cmd
    status = os.system(cmd)
    return status

if __name__ == '__main__':
    import os
    #result = pool.apply_async(f, [10])    # evaluate "f(10)" asynchronously
    #print result.get(timeout=2)           # prints "100" unless your computer is *very* slow
    exe = './mo_analyse'
    dir_target = ''
    sub_dir_s = [\
      #'HME_0.01/15.0_1.0E+04_1.0E-07_15.0_1.0E-02/',
      #'HME_0.01/30.0_1.0E+04_1.0E-07_15.0_1.0E-02/',
      #'HME_0.01/30.0_1.0E+06_1.0E-07_15.0_1.0E-02/',
      #'HME_0.01/15.0_1.0E+06_1.0E-07_15.0_1.0E-02/',
      'HME_0.01/15.0_1.0E+05_1.0E-07_15.0_1.0E-02/',
      'HME_0.01/30.0_1.0E+05_1.0E-07_15.0_1.0E-02/',
      'HME_0.01/20.0_1.0E+06_1.0E-07_15.0_1.0E-02/',
      'HME_0.01/20.0_1.0E+05_1.0E-07_15.0_1.0E-02/',
      'HME_0.01/20.0_1.0E+04_1.0E-07_15.0_1.0E-02/',
      'HME_0.01/10.0_1.0E+06_1.0E-07_15.0_1.0E-02/',
      'HME_0.01/10.0_1.0E+05_1.0E-07_15.0_1.0E-02/',
      #'HME_0.01/10.0_1.0E+04_1.0E-07_15.0_1.0E-02/',
      'HME_0.01/10.0_1.0E+03_1.0E-07_15.0_1.0E-02/',
      'HME_0.01/20.0_1.0E+03_1.0E-07_15.0_1.0E-02/',
      'HME_0.01/15.0_1.0E+03_1.0E-07_15.0_1.0E-02/',
      'HME_0.01/30.0_1.0E+03_1.0E-07_15.0_1.0E-02/',
      ]
    config_file = 'config_moment_.dat'
    species_to_analyse = 'SpeciesToAna.dat'
    cmd_s = [exe + ' ' + os.path.join(dir_target, subd, config_file) + ' ' + species_to_analyse for subd in sub_dir_s]

    pool = Pool(processes=2)
    print pool.map(f, cmd_s)


