# SCRaMbLE simulation #3 for circular chromosomes
import random
from comparison_sol import LoxP_unit_count_list

# This function find what is the smaller SV in a circular chromosome
def smallest_SV(pos1, pos2, syn_chr):
    LEN = len(syn_chr)
    # I want that the position 1 is always smaller than position 2
    if pos1 > pos2:
        temp = pos1
        pos1 = pos2
        pos2 = temp
    SV_fw = pos2 - pos1
    SV_rv = LEN - pos2 + 0 + pos1
    if SV_fw <= SV_rv:
        return True
    else:
        return False

#print(smallest_SV(2, 5, syn_chr))
#print(smallest_SV(5, 5, syn_chr))
#print(smallest_SV(5, 35, syn_chr))
#print("---------")

# If some essential LUs are not in the chromosome before SCRaMbLE, it could give problems.
def deletion_essential_circular(pos1, pos2, syn_chr, essential=[]):
    if pos1 > pos2:
        temp=pos1
        pos1=pos2
        pos2=temp
    # It check if all the essential fragments are in the synthetic chromosome before the deletion happen. If there are not, it simply continue with the deletion.
    syn_chr_abs = [abs(ele) for ele in syn_chr]
    T1 = all(elem in syn_chr_abs for elem in essential)

    # It create the new chromosome with the deletion
    if smallest_SV(pos1, pos2, syn_chr):
        new_chr = syn_chr[:pos1] + syn_chr[pos2:]
    else:
        new_chr = syn_chr[pos1:pos2]
    # It check if all the essential fragments are in the synthetic chromosome after the deletion happen. If there are not, it will output the original chr.
    new_chr_abs = [abs(ele) for ele in new_chr]
    T2 = all(elem in new_chr_abs for elem in essential)
    if T2 == False and T1 == True and essential != []:
        return syn_chr
    else:
        return new_chr

#print(deletion_essential_circular(3,7,syn_chr))
#print(deletion_essential_circular(2,1,syn_chr))
#print(deletion_essential_circular(5,35,syn_chr))

def inversion_circular(pos1, pos2, syn_chr):
    if pos1 > pos2:
        temp=pos1
        pos1=pos2
        pos2=temp
    if smallest_SV(pos1, pos2, syn_chr):
        inv_seq1 = syn_chr[pos1:pos2]
        inv_seq2 = []
        for loxP_f in inv_seq1:
            inv_seq2.insert(0, -loxP_f)
        return syn_chr[:pos1] + inv_seq2 + syn_chr[pos2:]
    else:
        inv_seq1 = syn_chr[pos2:] + syn_chr[:pos1]
        inv_seq2 = []
        for loxP_f in inv_seq1:
            inv_seq2.insert(0, -loxP_f)
        return inv_seq2[pos1:] + syn_chr[pos1:pos2] + inv_seq2[:pos1]

#print(inversion_circular(5,2,syn_chr))
#print(inversion_circular(35,5,syn_chr))
#print(inversion_circular(5,35,syn_chr))


def duplication_circular(pos1, pos2, syn_chr, CEN=[]):
    if pos1 > pos2:
        temp = pos1
        pos1 = pos2
        pos2 = temp
    if smallest_SV(pos1, pos2, syn_chr):        # It duplicates always the smaller portion of the chromosome
        new_chr = syn_chr[:pos1] + 2*syn_chr[pos1:pos2] + syn_chr[pos2:]
        if LoxP_unit_count_list(new_chr, CEN) > 1:  # The centromere is duplicated
            return syn_chr
        else:
            return new_chr
    else:
        duplicated = syn_chr[pos2:] + syn_chr[:pos1]
        new_chr = syn_chr[:pos1] + duplicated + syn_chr[pos1:]
        if LoxP_unit_count_list(new_chr, CEN) > 1:  # The centromere is duplicated
            return syn_chr
        else:
            return new_chr



def SCRaMbLE4_circular(syn_chr, Number_events, essential=[], mu=0, sigma=10, CEN=[], probability=[3, 2, 2, 1]):
    if len(syn_chr) < 2:
        return syn_chr
    # Add the Centromere to the essential LU
    for LU in CEN:
        if LU not in essential:
            essential.append(LU)
    new_chr = syn_chr
    events = random.choices(["NULL", "DEL", "INV", "DUP"], probability, k=Number_events)
    for event in events:
        pos1 = random.randrange(1, len(new_chr))
        # pick a random number to decide the length of the fragment. It is a Gaussian distribution, mu is the mean and sigma is the standard deviation.
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
        if random.random() >= 0.5:
            pos2 = pos1 + R
        else:
            pos2 = pos1 - R

        if pos2 >= len(new_chr):    # E.g. pos2=12 len(syn_chr)=10, new_pos2=12-10=2
            pos2 = pos2 % len(new_chr)
        if pos2 < 0:                # E.g. pos2=-2 len(syn_chr)=10, new_pos2=10-2=8
            if abs(pos2) > len(new_chr):    # E.g. pos2=-12 len(syn_chr)=10, new_pos2=-2, new_pos2=10-2=8
                pos2 = -(abs(pos2) % len(new_chr))
            pos2 = len(new_chr) + pos2
        if pos1 == pos2:            # Nothing happens
            continue
        # I want that the position 1 is always smaller than position 2
        if pos1 > pos2:
            temp = pos1
            pos1 = pos2
            pos2 = temp

        if pos1 < 0 or pos2 < 0 or pos1 > len(new_chr) or pos2 > len(new_chr):
            print("There is an error with the SCRaMbLE function. pos1 or pos2 are out of range.")
            print(new_chr)
            print("chr length =", len(new_chr))
            print("pos1 = ", pos1, "pos2 = ", pos2)

        if event == "NULL":
            new_chr = new_chr
        if event == "DEL":
            new_chr = deletion_essential_circular(pos1, pos2, new_chr, essential=essential)
        if event == "INV":
            new_chr = inversion_circular(pos1, pos2, new_chr)
        if event == "DUP":
            new_chr = duplication_circular(pos1, pos2, new_chr, CEN=CEN)
        # print(event, "between LoxP:", pos1, "and LoxP:", pos2)
        # print(new_chr)
    return new_chr

# test the code
if __name__ == "__main__":
    segments = 44  #number of loxP segments
    syn_chr = list(range(1, segments+1, 1))
    essential = [2,7,9,10,12,20]
    print("syn_chr =", syn_chr)
    for i in range(10):
        print(SCRaMbLE4_circular(syn_chr, 1000, essential=essential, CEN=[2]))
        print()
