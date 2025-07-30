import pandas as pd
import numpy as np
import time
import copy
import matplotlib.pyplot as plt
import plotly.figure_factory as ff
import datetime

# Load processing time data
ptd = pd.read_excel('datasheet.xlsx', sheet_name="Processing Time")
mcs = pd.read_excel('datasheet.xlsx', sheet_name="Machines Sequence")

# Convert data into list of lists
pt = []
mc = []
for i in range(len(ptd)):
    pt.append([None if pd.isna(x) else int(x) for x in ptd.iloc[i][1:]])
    mc.append([None if pd.isna(x) else int(x) for x in mcs.iloc[i][1:]])

# Calculate number of jobs and max number of machines in any job
num_jobs = len(pt)
num_mch = max(len([x for x in sq if x is not None]) for sq in mc)
nm_gene = sum(len([x for x in sq if x is not None]) for sq in mc)

# Input from user
population_size = int(input('Enter the size of population: ') or 20)
crss_rate = float(input('Enter the size of cross rate: ') or 0.9)
mut_rate = float(input('Enter the size of mutation rate: ') or 0.1)
mut_sel_rate = float(input('Enter the size of mutation selection rate: ') or 0.1)
nm_mut_jobs = round(nm_gene * mut_sel_rate)
nm_itr = int(input('Enter the number of iterations: ') or 1000)

start = time.time()

# Generate initial population
Tbest = float('inf')
bst_lst, bst_obj = [], []
population_lst = []
makespan_rec = []

for i in range(population_size):
    nm_random = list(np.random.permutation(nm_gene))
    population_lst.append(nm_random)
    for j in range(nm_gene):
        population_lst[i][j] = population_lst[i][j] % num_jobs

for m in range(nm_itr):
    Tbest_nw = float('inf')

    # 2-point crossover
    parent_lst = copy.deepcopy(population_lst)
    offspring_lst = copy.deepcopy(population_lst)
    s = list(np.random.permutation(population_size))

    for n in range(int(population_size / 2)):
        crs_prob = np.random.rand()
        if crss_rate >= crs_prob:
            p1 = population_lst[s[2 * n]][:]
            p2 = population_lst[s[2 * n + 1]][:]
            ch1 = p1[:]
            ch2 = p2[:]
            cut_point = list(np.random.choice(nm_gene, 2, replace=False))
            cut_point.sort()
            ch1[cut_point[0]:cut_point[1]] = p2[cut_point[0]:cut_point[1]]
            ch2[cut_point[0]:cut_point[1]] = p1[cut_point[0]:cut_point[1]]
            offspring_lst[s[2 * n]] = ch1[:]
            offspring_lst[s[2 * n + 1]] = ch2[:]

    # Repairment
    for m in range(population_size):
        job_count = {}
        larger, less = [], []

        for i in range(num_jobs):
            if i in offspring_lst[m]:
                count = offspring_lst[m].count(i)
                pos = offspring_lst[m].index(i)
                job_count[i] = [count, pos]
            else:
                count = 0
                job_count[i] = [count, 0]

            if count > len([x for x in mc[i] if x is not None]):
                larger.append(i)
            elif count < len([x for x in mc[i] if x is not None]):
                less.append(i)

        for a in range(len(larger)):
            chg_job = larger[a]
            while job_count[chg_job][0] > len([x for x in mc[chg_job] if x is not None]):
                for d in range(len(less)):
                    if job_count[less[d]][0] < len([x for x in mc[less[d]] if x is not None]):
                        offspring_lst[m][job_count[chg_job][1]] = less[d]
                        job_count[chg_job][1] = offspring_lst[m].index(chg_job)
                        job_count[chg_job][0] -= 1
                        job_count[less[d]][0] += 1
                    if job_count[chg_job][0] == len([x for x in mc[chg_job] if x is not None]):
                        break

    # Mutation
    for m in range(len(offspring_lst)):
        mut_prob = np.random.rand()
        if mut_rate >= mut_prob:
            m_chg = list(np.random.choice(nm_gene, nm_mut_jobs, replace=False))
            val_last = offspring_lst[m][m_chg[0]]
            for i in range(nm_mut_jobs - 1):
                offspring_lst[m][m_chg[i]] = offspring_lst[m][m_chg[i + 1]]
            offspring_lst[m][m_chg[nm_mut_jobs - 1]] = val_last

    # Calculate fitness value
    tot_chroms = copy.deepcopy(parent_lst) + copy.deepcopy(offspring_lst)
    ch_fitness, chroms_fit = [], []
    tot_fit = 0

    for m in range(population_size * 2):
        j_keys = [j for j in range(num_jobs)]
        key_count = {key: 0 for key in j_keys}
        j_count = {key: 0 for key in j_keys}

        m_keys = [j + 1 for j in range(num_mch)]
        m_availability = {key: 0 for key in m_keys}

        for i in tot_chroms[m]:
            if key_count[i] < len(mc[i]) and mc[i][key_count[i]] is not None:
                gen_t = int(pt[i][key_count[i]])
                gen_m = int(mc[i][key_count[i]])
                start_time = max(j_count[i], m_availability[gen_m])
                end_time = start_time + gen_t
                j_count[i] = end_time
                m_availability[gen_m] = end_time
                key_count[i] += 1

        makespan = max(j_count.values())
        ch_fitness.append(1 / makespan)
        chroms_fit.append(makespan)
        tot_fit += ch_fitness[m]

    # Selection (roulette wheel)
    pk, qk = [], []
    for i in range(population_size * 2):
        pk.append(ch_fitness[i] / tot_fit)
    for i in range(population_size * 2):
        cumulative = 0
        for j in range(0, i + 1):
            cumulative += pk[j]
        qk.append(cumulative)
    sel_rand = [np.random.rand() for i in range(population_size)]
    for i in range(population_size):
        if sel_rand[i] <= qk[0]:
            population_lst[i] = copy.deepcopy(tot_chroms[0])
        else:
            for j in range(0, population_size * 2 - 1):
                if sel_rand[i] > qk[j] and sel_rand[i] <= qk[j + 1]:
                    population_lst[i] = copy.deepcopy(tot_chroms[j + 1])
                    break

    # Comparison
    for i in range(population_size * 2):
        if chroms_fit[i] < Tbest_nw:
            Tbest_nw = chroms_fit[i]
            sq_nw = copy.deepcopy(tot_chroms[i])
        if Tbest_nw <= Tbest:
            Tbest = Tbest_nw
            sq_best = copy.deepcopy(sq_nw)

    makespan_rec.append(Tbest)

# Output results
print("Optimal sequence", sq_best)
print("Optimal value: %f" % Tbest)
print('Elapsed time: %s' % (time.time() - start))

# Plot makespan over generations
plt.plot([i for i in range(len(makespan_rec))], makespan_rec, 'b')
plt.ylabel('Makespan', fontsize=15)
plt.xlabel('Generation', fontsize=15)
plt.show()

# Plot Gantt chart
m_keys = [j + 1 for j in range(num_mch)]
j_keys = [j for j in range(num_jobs)]
key_count = {key: 0 for key in j_keys}
j_count = {key: 0 for key in j_keys}
m_count = {key: 0 for key in m_keys}
j_record = {}

for i in sq_best:
    if key_count[i] < len(mc[i]) and mc[i][key_count[i]] is not None:
        gen_t = int(pt[i][key_count[i]])
        gen_m = int(mc[i][key_count[i]])
        start_time = max(j_count[i], m_count[gen_m])
        end_time = start_time + gen_t
        j_count[i] = end_time
        m_count[gen_m] = end_time

        start_time_str = str(datetime.timedelta(seconds=start_time))
        end_time_str = str(datetime.timedelta(seconds=end_time))

        j_record[(i, gen_m)] = [start_time_str, end_time_str]

        key_count[i] += 1

df = []
for m in m_keys:
    for j in j_keys:
        if (j, m) in j_record:
            df.append(dict(Task='Machine %s' % (m), Start='2024-05-17 %s' % (str(j_record[(j, m)][0])),
                           Finish='2024-05-17 %s' % (str(j_record[(j, m)][1])), Resource='Job %s' % (j + 1)))

print("Gantt chart data:")
print(df)

# Create Gantt chart
fig = ff.create_gantt(df, index_col='Resource', show_colorbar=True, group_tasks=True, showgrid_x=True, title='Job Shop Schedule')

# Write Gantt chart to HTML file and open it in default web browser
fig.write_html('gantt_chart.html', auto_open=True)
