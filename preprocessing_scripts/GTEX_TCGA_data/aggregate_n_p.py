import csv, glob

base_dir = '/home/CAM/adamson/Recount2_files/'

#regression_outputs_trad
regression_outputs_trad = glob.glob(base_dir + 'regression_outputs_trad/*_predictions.tsv')
summary_stats = {}
d = open(base_dir + 'n_p_reg_summary.tsv', 'w')
writer = csv.writer(d, delimiter = '\t')
writer.writerow(['reg_id', 'n', 'p_dummy', 'p_trad_genes', 'p_select_genes']) 
for file_ in regression_outputs_trad:
    #sample_id  observed_PSI    predicted_PSI
    #50099   0.9478599999999999  0.9572223206658366
    reg_id = file_.replace(base_dir, '').replace('regression_outputs_trad/', '').replace('_predictions.tsv', '')
    a = open(file_, 'r')
    n = -1 #header
    for line in a:
        n +=1
    a.close()

    b = open(file_.replace('_predictions', '_coefficients'), 'r')
    reader = csv.DictReader(b, delimiter = '\t')
    #coefficient    std_error   t_val   p_val   variable_id
    #0.242914097046428   3.317473762878164   0.07322261287024658 0.9416299690640935  (Intercept)
    #-0.7881439949122776 0.37417300518551727 -2.106362522121314  0.03518979002451379 batch10
    p_dummy = 0; p_trad = 0
    for line in reader:
        if 'ENSG' not in line['variable_id'] and line['variable_id'] != '(Intercept)':
            p_dummy += 1
        elif 'ENSG' in line['variable_id']:
            p_trad += 1
    b.close()

    c = open(file_.replace('_predictions', '_coefficients').replace('_trad', ''), 'r')
    reader = csv.reader(c, delimiter = '\t')
    p_select = 0
    for line in reader:
        if 'ENSG' in line[0]:
            p_select += 1
    c.close()
    writer.writerow([reg_id, n, p_dummy, p_trad, p_select])
d.close()

