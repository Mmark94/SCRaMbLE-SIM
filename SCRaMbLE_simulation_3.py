# SCRaMbLE simulation 3
import random
import matplotlib.pyplot as plt
from SCRaMbLE_simulation_3_circular import SCRaMbLE4_circular
#from Mapping_coverage_MM import plot_LU_CN
#from Mapping_coverage_MM import plot_LU_CN_percentage

# This function was copied from comparison_sol.py
def LoxP_unit_count_list(Path, list_unit):
    LP_unit_count = 0
    for unit in list_unit:
        LP_unit_count = LP_unit_count + Path.count(unit) + Path.count(-unit)
    return LP_unit_count

# This function was copied from comparison_sol.py
def invert(x):
    inv_x=[]
    for element in x:
        if isinstance(element, int):
            inv_x.insert(0,-element)
        else:
            inv_x.insert(0,element)
    return inv_x

def deletion(pos1, pos2, syn_chr):
    if pos1 > pos2:
        temp=pos1
        pos1=pos2
        pos2=temp
    return syn_chr[:pos1] + syn_chr[pos2:]

def deletion_essential0(pos1, pos2, syn_chr, essential=[]):
    if pos1 > pos2:
        temp=pos1
        pos1=pos2
        pos2=temp
    # It check if all the essential fragments are in the synthetic chromosome before the deletion happen. If there are not, it simply continue with the deletion.
    syn_chr_abs = [abs(ele) for ele in syn_chr]
    T1 = all(elem in syn_chr_abs for elem in essential)
    # It create the new chromosome with the deletion
    new_chr = syn_chr[:pos1] + syn_chr[pos2:]
    # It check if all the essential fragments are in the synthetic chromosome after the deletion happen. If there are not, it will output the original chr.
    new_chr_abs = [abs(ele) for ele in new_chr]
    T2 = all(elem in new_chr_abs for elem in essential)
    if T2 is False and T1 is True and essential != []:
        return syn_chr
    else:
        return new_chr

def deletion_essential(pos1, pos2, syn_chr, essential=[]):
    if pos1 > pos2:
        temp=pos1
        pos1=pos2
        pos2=temp
    # It create the new chromosome with the deletion
    new_chr = syn_chr[:pos1] + syn_chr[pos2:]
    if essential==[]:
        return new_chr
    # It check if all the essential fragments are in the synthetic chromosome before the deletion happen. If there are not, it simply continue with the deletion.
    new_essential = []
    syn_chr_abs = [abs(ele) for ele in syn_chr]
    for esse in essential:
        if esse in syn_chr_abs:
            new_essential.append(esse)
    # It check if all the essential fragments are in the synthetic chromosome after the deletion happened. If there are not, it will output the original chr.
    new_chr_abs = [abs(ele) for ele in new_chr]
    T = all(elem in new_chr_abs for elem in new_essential)
    if T is False:
        return syn_chr
    else:
        return new_chr

# Synthetic lethal interactions. These are genes that are redundant, they are non-essentials but if another genes is missing, they could become essentials.

def syn_lethal_interactions(syn_chr, pairs=[[]]):
    for P in pairs:
        # If at least one of the LUs in the pair P is present in the chromosome syn_chr, so the cell is viable
        viable = not(all(elem not in syn_chr for elem in P))
        print(P, viable)
        if viable is False:
            return False
    return True
#syn_chr = [1,2,3,4,5,6,7]
#pairs = [[2,5], [3,1,5], [4,8], [8,10]]
#print(syn_lethal_interactions(syn_chr, pairs=pairs))

def inversion(pos1, pos2, syn_chr):
    if pos1 > pos2:
        temp=pos1
        pos1=pos2
        pos2=temp
    return syn_chr[:pos1] + invert(syn_chr[pos1:pos2]) + syn_chr[pos2:]


def duplication(pos1, pos2, syn_chr, CEN=[]):
    if pos1 > pos2:
        temp=pos1
        pos1=pos2
        pos2=temp
    new_chr = syn_chr[:pos1] + 2*syn_chr[pos1:pos2] + syn_chr[pos2:]
    if LoxP_unit_count_list(new_chr, CEN) > 1:
        return syn_chr
    else:
        return new_chr

# When duplications happen they could be tandem duplication or inverted duplication.
def inverted_duplication(pos1, pos2, syn_chr, CEN=[]):
    if pos1 > pos2:
        temp=pos1
        pos1=pos2
        pos2=temp
    new_chr = syn_chr[:pos1] + syn_chr[pos1:pos2] + invert(syn_chr[pos1:pos2]) + syn_chr[pos2:]
    if LoxP_unit_count_list(new_chr, CEN) > 1:
        return syn_chr
    else:
        return new_chr

#A = [1,2,3,4,5,6]
#print(deletion_essential(0,3,A, [3,4,7]))
#print(deletion(3,0,A))
#print(inversion(3,0,A))
#print(inverted_duplication(3, 0, A))

# This is a reciprocal translocation. This means that no LUs are lost.
# However, it might happen that the centromere is translocated. In this case, there might be a chromosome with 0 or 2 centromeres (not viable).
def translocation(chr1, chr2):
    if len(chr1) <= 1 or len(chr2) <= 1:
        return [chr1, chr2]
    pos1 = random.randrange(1, len(chr1))
    pos2 = random.randrange(1, len(chr2))
    new_chr1 = chr1[:pos1] + chr2[pos2:]
    new_chr2 = chr2[:pos2] + chr1[pos1:]
    return [new_chr1, new_chr2]

#A = [1,2,3,4,5,6]
#B = [12,13,14,15,16,17,18,19,20]
#print(translocation(A , B))

def SCRaMbLE(syn_chr):
    pos1 = random.randrange(1, len(syn_chr))
    pos2 = random.randrange(1, len(syn_chr))
    if pos1 > pos2:
        temp=pos1
        pos1=pos2
        pos2=temp
    if pos1 == pos2:
        print("Nothing happened")
        return syn_chr
    event = random.choices(["NULL", "DEL", "INV", "DUP"], [3, 2, 2, 1], k=1)[0]
    if event == "NULL":
        return syn_chr
    print(event, "between LoxP:", pos1, "and LoxP:", pos2)
    if event == "DEL":
        return deletion(pos1, pos2, syn_chr)
    if event == "INV":
        return inversion(pos1, pos2, syn_chr)
    if event == "DUP":
        return duplication(pos1, pos2, syn_chr)



def repeated_SCRaMbLE(syn_chr, Number_events):
    new_chr = syn_chr
    for event in range(Number_events):
        new_chr = SCRaMbLE(new_chr)
        print(new_chr)
    return new_chr



def SCRaMbLE2(syn_chr, Number_events):
    new_chr = syn_chr
    events = random.choices(["NULL", "DEL", "INV", "DUP"], [3, 2, 2, 1], k=Number_events)
    for event in events:
        pos1 = random.randrange(2, len(syn_chr))
        pos2 = random.randrange(2, len(syn_chr))
        if pos1 > pos2:
            temp = pos1
            pos1 = pos2
            pos2 = temp
        if pos1 == pos2:
            #print("Nothing happened")
            continue
        if event == "NULL":
            new_chr = new_chr
        if event == "DEL":
            new_chr = deletion(pos1, pos2, new_chr)
        if event == "INV":
            new_chr = inversion(pos1, pos2, new_chr)
        if event == "DUP":
            new_chr = duplication(pos1, pos2, new_chr)
        #print(event, "between LoxP:", pos1, "and LoxP:", pos2)
        #print(new_chr)
    return new_chr


def SCRaMbLE3(syn_chr, Number_events, mu=7, sigma=7):
    new_chr = syn_chr
    events = random.choices(["NULL", "DEL", "INV", "DUP"], [3, 2, 2, 1], k=Number_events)
    for event in events:
        pos1 = random.randrange(2, len(syn_chr))
        # pick a random number to decide the length of the fragment
        R = int(random.gauss(mu, sigma))
        if R < 1:
            R = 1
        pos2 = pos1 + R
        if pos2 > len(syn_chr):
            pos2 = len(syn_chr)
        if event == "NULL":
            new_chr = new_chr
        if event == "DEL":
            new_chr = deletion(pos1, pos2, new_chr)
        if event == "INV":
            new_chr = inversion(pos1, pos2, new_chr)
        if event == "DUP":
            new_chr = duplication(pos1, pos2, new_chr)
        #print(event, "between LoxP:", pos1, "and LoxP:", pos2)
        #print(new_chr)
    return new_chr

def SCRaMbLE4(syn_chr, Number_events, essential=[], mu=0, sigma=10, CEN=[], probability=[3, 2, 2, 1], event_type=False):
    if len(syn_chr) == 0:
        if event_type:
            return syn_chr, []
        return syn_chr
    # Add the Centromere to the essential LU
    for LU in CEN:
        if LU not in essential:
            essential.append(LU)

    new_chr = syn_chr[:]
    events = random.choices(["NULL", "DEL", "INV", "DUP"], probability, k=Number_events)
    for event in events:
        pos1 = random.randrange(0, len(new_chr), 1)
        # pick a random number to decide the length of the fragment. In the Gaussian distribution, mu is the mean and sigma is the standard deviation.
        # Discretized truncated normal distribution (DTND)
        #R = int(random.gauss(mu, sigma))
        #if R < 1:
        #    R = 1
        # Half-normal distribution. https://en.wikipedia.org/wiki/Half-normal_distribution
        R = 0
        while R == 0:
            # this loop makes repeat the choosing if R == 0
            R = abs(int(random.gauss(mu, sigma)))

        # the second cut position pos2 could be before or after the first position pos1
        if random.random() >= 0.5 or pos1 == 0:
            pos2 = pos1 + R
        else:
            pos2 = pos1 - R

        if pos2 >= len(new_chr):
            pos2 = len(new_chr)
        if pos2 < 0:
            pos2 = 0
        if pos1 == pos2:    # Nothing happens
            continue
        # I want that the position 1 is always smaller than position 2
        if pos1 > pos2:
            temp = pos1
            pos1 = pos2
            pos2 = temp

        # This is just to catch some errors in the function
        if pos1 < 0 or pos2 < 0 or pos1 > len(new_chr) or pos2 > len(new_chr):
            print("There is an error with the SCRaMbLE function. pos1 or pos2 are out of range.")
            print(new_chr)
            print("chr length =", len(new_chr))
            print("pos1 = ", pos1, "pos2 = ", pos2)

        if event == "NULL":
            continue
        if event == "DEL":
            #new_chr = deletion(pos1, pos2, new_chr)
            new_chr = deletion_essential(pos1, pos2, new_chr, essential)
        if event == "INV":
            new_chr = inversion(pos1, pos2, new_chr)
        if event == "DUP":
            new_chr = duplication(pos1, pos2, new_chr, CEN)
        #print(event, "between LoxP:", pos1, "and LoxP:", pos2)
        #print(new_chr)
    if event_type:
        return new_chr, events
    return new_chr

# Plot the SCRaMbLE events length. This should follow a "discretized half normal distribution"
def plot_events_length(events=100000, mu=0, sigma=10):
    syn_chr = list(range(1, 45))
    Dict_E_length = {}
    for i in range(len(syn_chr)+1):
        Dict_E_length[i] = 0

    for _ in range(events):
        pos1 = random.randrange(0, len(syn_chr), 1)
        # pick a random number to decide the length of the fragment. In the Gaussian distribution, mu is the mean and sigma is the standard deviation.
        # Discretized truncated normal distribution (DTND)
        #R = int(random.gauss(mu, sigma))
        #if R < 1:
        #    R = 1
        # Half-normal distribution
        R = 0
        while R == 0:
            R = abs(int(random.gauss(mu, sigma)))

        # the second cut position pos2 could be before or after the first position pos1
        if random.random() >= 0.5 or pos1 == 0:
            pos2 = pos1 + R
        else:
            pos2 = pos1 - R

        if pos2 >= len(syn_chr):
            pos2 = len(syn_chr)
        if pos2 < 0:
            pos2 = 0
        E_length = abs(pos2 - pos1)
        Dict_E_length[E_length] = Dict_E_length[E_length] + 1
    # remove keys with zero events.
    for i in range(len(syn_chr)+1):
        if Dict_E_length[i] == 0:
            Dict_E_length.pop(i, None)
    #Dict_E_length.pop(1, None)
    #Dict_E_length[1] = Dict_E_length[2] + 150

    # Convert the event length number in probability
    Dict_E_length_probability = {}
    for Key, Value in Dict_E_length.items():
        Dict_E_length_probability[Key] = Value / events
    print("Dict_E_length_probability =", Dict_E_length_probability)
    cumulative_p = 0
    N50 = 0
    N75 = 0
    for Key, Value in Dict_E_length_probability.items():
        cumulative_p = cumulative_p + Value
        print(Key, Value, cumulative_p)
        if cumulative_p > 0.5 and N50 == 0:
            N50 = Key - 0.5
        if cumulative_p > 0.75 and N75 == 0:
            #N75 = Key - 0.5
            N75 = Key + 0.5
    # Plot the Events length
    plt.figure(figsize=(7, 3.5), dpi=300)
    # Font size
    SMALL_SIZE = 8
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 11
    plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
    plt.grid(True, axis="y", zorder=-1, alpha=0.5)
    plt.bar(Dict_E_length.keys(), Dict_E_length_probability.values(), align='center', zorder=2)
    plt.vlines(N50, 0, max(Dict_E_length_probability.values()) * 0.95, colors="orange", linestyle=":", label="N50", zorder=1)
    plt.vlines(N75, 0, max(Dict_E_length_probability.values()) * 0.95, colors="green", linestyle=":", label="N75", zorder=1)
    plt.xticks(range(max(Dict_E_length.keys()) + 1), range(max(Dict_E_length_probability.keys()) + 1))
    plt.ylabel("Event length probability")
    plt.xlabel("Length (LUs)")
    plt.title("SCRaMbLE events length distribution", fontsize=10)
    #plt.text(max(Dict_E_length.keys()) * 0.65, max(Dict_E_length_probability.values()) * 0.78, "Number of Events = " + str(events) + "\n" + "Mean length (mu) = " + str(mu) + "\n" + "Sigma = " + str(sigma) + " LUs")
    plt.text(20, max(Dict_E_length_probability.values()) * 0.78, "Number of Events = " + str(events) + "\n" + "Mean length (mu) = " + str(mu) + "\n" + "Sigma = " + str(sigma) + " LUs")
    plt.xlim((0, 32))
    #plt.xticks(rotation=90)
    plt.savefig("Events_length_distribution.png", dpi=300, bbox_inches='tight')
    plt.savefig("Events_length_distribution.svg", format='svg', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close()
    return Dict_E_length_probability

# This make sure that there is at last one SCRaMbLE event each time
def force_SCRaMLE(syn_chr, Number_events, essential=[], mu=0, sigma=7, CEN=[], probability=[0, 2, 2, 1]):
    new_chr2 = syn_chr[:]
    for _ in range(Number_events):
        new_chr1 = new_chr2[:]
        counter = 0         # the counter make sure the program do not get stuck and loop infinite times.
        while new_chr1 == new_chr2 and counter < 20:
            new_chr2 = SCRaMbLE4(new_chr1, 1, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability)
            counter = counter + 1
        #print(new_chr2)
    return new_chr2

# This make sure that there is at last one SCRaMbLE event each time
def force_SCRaMLE_lin_cir(syn_chr: list, Number_events: int, essential=[], circular=False, mu=0, sigma=7, CEN=[], force=True, probability=[0, 2, 2, 1]):
    if not(force):
        if circular:
            return SCRaMbLE4_circular(syn_chr=syn_chr, Number_events=Number_events, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability)
        else:
            return SCRaMbLE4(syn_chr=syn_chr, Number_events=Number_events, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability)
    new_chr2 = syn_chr[:]
    for _ in range(Number_events):
        new_chr1 = new_chr2[:]
        counter = 0         # the counter make sure the program do not get stuck and loop infinite times.
        while new_chr1 == new_chr2 and counter < 20:
            if circular:
                new_chr2 = SCRaMbLE4_circular(new_chr1, 1, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability)
            else:
                new_chr2 = SCRaMbLE4(new_chr1, 1, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability)
            counter = counter + 1
        #print(new_chr2)
    return new_chr2

# This records also the event type
def force_SCRaMLE_lin_cir_events(syn_chr: list, Number_events: int, essential=[], circular=False, mu=0, sigma=7, CEN=[], force=True, probability=[0, 2, 2, 1], event_type=False):
    if not(force):
        if circular:
            return SCRaMbLE4_circular(syn_chr=syn_chr, Number_events=Number_events, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability)
        else:
            return SCRaMbLE4(syn_chr=syn_chr, Number_events=Number_events, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability)
    new_chr2 = syn_chr[:]
    Events = []
    for _ in range(Number_events):
        new_chr1 = new_chr2[:]
        counter = 0         # the counter make sure the program do not get stuck and loop infinite times.
        while new_chr1 == new_chr2 and counter < 20 and new_chr2 != []:
            if event_type:
                if circular:
                    chr_temp = SCRaMbLE4_circular(new_chr1, 1, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability, event_type=True)
                else:
                    chr_temp = SCRaMbLE4(new_chr1, 1, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability, event_type=True)
                new_chr2 = chr_temp[0]
                if new_chr1 != new_chr2:
                    Events.append(chr_temp[1][0])
            else:
                if circular:
                    new_chr2 = SCRaMbLE4_circular(new_chr1, 1, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability)
                else:
                    new_chr2 = SCRaMbLE4(new_chr1, 1, essential=essential, mu=mu, sigma=sigma, CEN=CEN, probability=probability)
            counter = counter + 1
            if counter == 20:
                Events.append("NULL")
        #print(new_chr2)
    if event_type:
        #print("Events =", len(Events), Events)
        return new_chr2, Events
    return new_chr2

#syn_chr = list(range(1, 45, 1))
#essential = [2,7,9,10,12,20]
#print(force_SCRaMLE(syn_chr, 500, essential))

"""
mu=7
sigma=7
for i in range(50):
    R = int(random.gauss(mu, sigma))
    print(R)
"""

def SCRaMbLE4_lin_cir(syn_chr, Number_events, essential=[], circular=False, mu=7, sigma=7, CEN=[], probability=[3, 2, 2, 1]):
    if circular:
        return SCRaMbLE4_circular(syn_chr, Number_events, essential, mu, sigma, CEN=CEN, probability=probability)
    else:
        return SCRaMbLE4(syn_chr, Number_events, essential, mu, sigma, CEN=CEN, probability=probability)


# This function can SCRaMbLE multiple chromosomes. The input is a list of chromosomes. The probability of SCRaMbLE in one chromosome is proportional with its length and there is a small probabiity of translocation (Ptra=0.05).
def SCRaMbLE_muliple_chrs(list_chr: list, Number_events=1, essential=[], circular=False, mu=0, sigma=7, CEN=[], force=True, probability=[0, 2, 2, 1], Ptra=0.05):
    #Number_events = 1
    if isinstance(list_chr[0], list):       # There are multiple chromosomes
        new_chr = list_chr[:]
        # Name each chromosome starting from 0, 1, 2, ...
        list_chr_name = list(range(len(list_chr)))
        # Find the length of all the chromosomes
        list_chr_len = [len(x) for x in list_chr]
        # Decide where the SCRaMbLE event should happen. More the chr is long more is the likelihood to be chosen.
        chr_events = random.choices(list_chr_name, list_chr_len, k=Number_events)
        #print("list_chr_name =", list_chr_name)
        #print("list_chr_len =", list_chr_len)
        #print("chr_events =", chr_events)
        for i in range(Number_events):
            #print(i)
            # Decide if the SCRaMbLE event is intra_chr or inter_chr (translocation). The translocations have a probability of 5% to happen.
            if random.random() < Ptra:
                #print("Translocation")
                # Keep looping until the second translocated chromosome is different from the first. Maximum for 10 iterations.
                for _ in range(10):
                    second_tra_chr = random.choices(list_chr_name, list_chr_len, k=1)[0]
                    if second_tra_chr != chr_events[i]:
                        break
                if second_tra_chr == chr_events[i]:
                    print("Null SCRaMbLE event because of translocation. Could not chose the chromosomes.")
                    continue
                # Generate the translocation between the two chromosomes chosen
                translocated_chr = translocation(new_chr[chr_events[i]], new_chr[second_tra_chr])
                # If the centromere list is not provided, do nothing
                if CEN == []:
                    new_chr[chr_events[i]] = translocated_chr[0]
                    new_chr[second_tra_chr] = translocated_chr[1]
                else:
                    # Check if the translocated chromosomes have 1 and only 1 CEN. If they have 0 or 2 or more CEN, the program will discard the translocation.
                    CEN_chr1 = LoxP_unit_count_list(Path=translocated_chr[0], list_unit=CEN)
                    CEN_chr2 = LoxP_unit_count_list(Path=translocated_chr[1], list_unit=CEN)
                    if CEN_chr1 == 1 and CEN_chr2 == 1:
                        new_chr[chr_events[i]] = translocated_chr[0]
                        new_chr[second_tra_chr] = translocated_chr[1]
                    else:
                        print("Null SCRaMbLE event because of translocation. Number of centromeres for chromosome =", CEN_chr1, CEN_chr2)
            else:
                new_chr[chr_events[i]] = force_SCRaMLE_lin_cir(new_chr[chr_events[i]], Number_events, essential=essential, circular=circular, mu=mu, sigma=sigma, CEN=CEN, force=force, probability=probability)
        return new_chr
    else:
        return force_SCRaMLE_lin_cir(list_chr, Number_events, essential=essential, circular=circular, mu=mu, sigma=sigma, CEN=CEN, force=force, probability=probability)

def SCRaMbLE_muliple_chrs_events(list_chr: list, Number_events=10, essential=[], circular=False, mu=0, sigma=7, CEN=[], force=True, probability=[0, 2, 2, 1]):
    new_chr = list_chr[:]
    if not (force):
        for _ in range(Number_events):
            new_chr = SCRaMbLE_muliple_chrs(new_chr, essential=essential, circular=circular, mu=mu, sigma=sigma, CEN=CEN, force=force, probability=probability)
    else:
        new_chr2 = new_chr[:]
        for _ in range(Number_events):
            new_chr = new_chr2[:]
            counter = 0         # the counter make sure the program do not get stuck and loop infinite times.
            while new_chr == new_chr2 and counter < 20:
                new_chr2 = SCRaMbLE_muliple_chrs(new_chr, essential=essential, circular=circular, mu=mu, sigma=sigma, CEN=CEN, force=force, probability=probability)
                counter += 1
    return new_chr

#A = [[1,2,3,4,5,6],[7,8,9,10,11],[12,13,14,15,16,17,18,19,20]]
#print(A)
#print(SCRaMbLE_muliple_chrs(A))
#print(SCRaMbLE_muliple_chrs_events(A, 20, [3,10,13], 7, 7, [3,10,13]))


# test the code
if __name__ == "__main__":

    segments = 44  #number of loxP segments
    syn_chr = list(range(1, segments+1, 1))
    essential = [2, 7, 9, 10, 12, 19, 20, 24]  # LUs 19 and 24 are not essential but required for fast growth. Deletion of LU 6 can also generate some slow growth phenotype.
    print("syn_chr =", syn_chr)

    for i in range(0):
        print(SCRaMbLE4_lin_cir(syn_chr, 10, essential, False, probability=[2,1,1,2]))
        print(SCRaMbLE4_lin_cir(syn_chr, 10, essential, True, probability=[2,1,1,2]))
        print(force_SCRaMLE_lin_cir(syn_chr, 10, essential, False, probability=[2,1,1,2]))
        print(force_SCRaMLE_lin_cir_events(syn_chr, 10, essential, False, probability=[2, 1, 1, 2], event_type=True))
        print()

    #plot_events_length(events=1000000, mu=0, sigma=10)
