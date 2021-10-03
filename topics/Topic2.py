# ---
# jupyter:
#   jupytext:
#     formats: ipynb,md,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.8.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# + [markdown] slideshow={"slide_type": "slide"}
# # Lab 2 - How do we compare DNA sequences?
# ## Dynamic Programming
# Material and embedded lab.
#
# Motivation and some exercises are variations on those available in Bioinformatics Algorithms: An Active-Learning Approach by Phillip Compeau & Pavel Pevzner.

# + slideshow={"slide_type": "skip"}

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Announcements
# GitHub Classroom and Canvas, Assignment on Canvas, etc
#
# Project update - I hope everyone is rocking and rolling. Please chime in now if you would like :)
#
# Data frames/signatures/etc - you may immediately translate the pandas structures that I uses into your own Python generic data structures. I find this to make the problems harder, but this could make the problems easier for you. It is up to you :) You must satisfy the autograder of course so don't change the input or output types.
#
# Assignments for week
#
# ## Plans
# See headings on slide deck
#
#

# + [markdown] slideshow={"slide_type": "subslide"}
# Learning objectives for the week and lab:
# 1. Build our mental model of biology by incorporating comparison as a framework for understanding biology
# 2. Apply your dynamic programming experience from algorithms
# 3. Gain experience translating a biological problem into a problem we can solve via code

# + [markdown] slideshow={"slide_type": "slide"}
# # History and motivation
#
# Searching all new sequences against sequence databases is now the first order of business in genomics!
#
# Why is it so important to us?

# + [markdown] slideshow={"slide_type": "subslide"}
# Mohamed Marahiel conjectured that since A-domains of protins (adenylation domains) have the same function (i.e., adding an amino acid to the growing peptide), different A-domains should have similar parts. Each A-domain is about 500 amino acids long and is responsible for adding a single amino acid.

# + [markdown] slideshow={"slide_type": "subslide"}
# Taking 3 common A-domains and putting them down one after another, there are only three conserved columns (shown in red below) are common to the three sequences and have likely arisen by pure chance:
#
# <img src="http://bioinformaticsalgorithms.com/images/Alignment/A_domain_2.png" width=2000>

# + [markdown] slideshow={"slide_type": "fragment"}
# Well that stinks...

# + [markdown] slideshow={"slide_type": "subslide"}
# If we slide the second sequence only one amino acid to the right, adding a space symbol ("-") to the beginning of the sequence, then we find 11 conserved columns!
#
# <img src="http://bioinformaticsalgorithms.com/images/Alignment/A_domain_3.png" width=2000>

# + [markdown] slideshow={"slide_type": "subslide"}
# Adding a few more space symbols reveals 14 conserved columns:
#
# <img src="http://bioinformaticsalgorithms.com/images/Alignment/A_domain_4.png" width=2000>

# + [markdown] slideshow={"slide_type": "fragment"}
# Now that is getting better, and the number of conserved columns is going up!

# + [markdown] slideshow={"slide_type": "subslide"}
# It turns out that the red columns represent the conserved core shared by many A-domains. Now that Marahiel knew how to correctly align the A-domains, he hypothesized that some of the remaining variable columns should code for Asp, Orn, and Val. He discovered that the non-ribosomal code is defined by 8 amino acid-long non-ribosomal signatures, which are shown as purple columns below.
#
# <img src="http://bioinformaticsalgorithms.com/images/Alignment/A_domain_6.png" width=2000>

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Quick aside about the book chapter
#
# So the book has a very interesting introduction to this topic, but I would say that it requires more biology than currently within reach of this class in a single pass (i.e., one week). We will be layering in some of this biology this week and we'll build on it throughout the semester.

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Identification of homologous genes
# * Two genes are homologous if they share a common ancestor
# * Living organisms share a large number of genes descended from common ancestors 
# * Functionality is preserved but differences in sequence accumulate as they diverge from each other. 
# * These differences may be due to mutations that change a symbol (nucleotide or amino acid) for another or insertions / deletions, indels, which insert or delete a symbol in the corresponding sequence.

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Our life this lab/lecture
#
# <img src="https://raw.githubusercontent.com/gregcaporaso/An-Introduction-To-Applied-Bioinformatics/master/book/fundamentals/images/alignment.png">

# + [markdown] slideshow={"slide_type": "subslide"}
# ## An important distinction (and a biologist will correct you every time)
# * Homology is a dichotomous characteristic, i.e., given two genes are either homologous genes or not. It is binary.
#
# * However, given two sequences corresponding to two genes, can be said that there are different levels of similarity based on an alignment between them. 
# * Our key question is to determine whether a good alignment between two sequences is significant enough to consider that both genes are homologous. 
#     * This task is done through a hypothesis testing and the corresponding p-values are used to make a decision.

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Two different forms of homology. 
# Paralogs - When the origin of two homologous genes is due to a process of gene duplication within the same species
#
# Orthologous genes - origin is due to a speciation process resulting in homologous genes in these different species

# + [markdown] slideshow={"slide_type": "slide"}
# # Introduction to Sequence Alignment
# How do you align/match up ATGCTTA and TGCATTAA that have subtle similarities?
# <pre>
# ATGC-TTA-
# -TGCATTAA
# </pre>
#
# It's all a game. The goal is to maximize "points" which for us is matching nucleotides or amino acids.
#
# At each turn, you have two choices
# 1. You can remove the first symbol from both sequences and align them. You'll earn a point if they match.
# 2. You can remove the first symbol from either of the two sequences in which case you earn no points, but you may set yourself up to earn more points later. 

# + [markdown] slideshow={"slide_type": "subslide"}
# ### Greedy approach
# Let's say that we want to take a greedy approach to alignment. Meaning we will only consider the choice in front of us. Example:
#
# <pre>
# s1=AACCTTGG
# s2=ACACTGTGA
# </pre>
#
# For each move in the game: 
#     * if s1[0] == s2[0], then add s1[0] to the longest common subsequence. 
#     * else randomly choose to either remove s1[0], s2[0], or both s1[0] and s2[0]

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 1** Use a greedy approach to return suboptimal (or optimal) solutions to the longest common subsequence problem.
#
# Find a longest common subsequence of two strings
#
# Input: Two strings
#
# Output: A longest common subsequence of these strings

# + slideshow={"slide_type": "subslide"}
import random

def random_action(s1,s2):
    assert len(s1) > 0 and len(s2) > 0
    choices = [(s1[1:],s2[1:]),(s1[1:],s2), (s1,s2[1:])]
    return random.choice(choices)
    
def greedy_lcs(s1,s2,seed=0):
    random.seed(seed)
    lcs = ""
    # YOUR SOLUTION HERE
    return lcs


print(greedy_lcs("AACCTTGG","ACACTGTGA",seed=0))
print(greedy_lcs("AACCTTGG","ACACTGTGA",seed=100))
print(greedy_lcs("AACCTTGG","ACACTGTGA",seed=1000))
print(greedy_lcs("AACCTTGG","ACACTGTGA",seed=2000))


# + [markdown] slideshow={"slide_type": "subslide"}
# Well... That was easy to implement, but the longest common subsequence is AACTGG, so we did not really solve the problem. Before we move on though, let's modify our algorithm so it returns the alignment (i.e., with indels and mutations).
#
# **Exercise 2** Modify your solution to exercise 1 to return the alignment as two strings with "-" characters when there is a indel. Do not modify ``random_action`` as that is what the autograder is going to rely on you using. HINT: Just keep track of the strings before calling ``random_action``.

# + slideshow={"slide_type": "subslide"}
def greedy_alignment(s1,s2,seed=0):
    random.seed(seed)
    s1_new = ""
    s2_new = ""
    return "\n".join([s1_new,s2_new])


print(greedy_alignment("AACCTTGG","ACACTGTGA",seed=0))
print()
print(greedy_alignment("AACCTTGG","ACACTGTGA",seed=100))
print()
print(greedy_alignment("AACCTTGG","ACACTGTGA",seed=1000))
print()
print(greedy_alignment("AACCTTGG","ACACTGTGA",seed=2000))

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Why do we have mismatches in our alignment? Why not just gaps?
# When consider both nucleotide sequences and amino acid sequences (different alphabets if you are a pure CS thinker), then sometimes you can achieve a better overall alignment if you allow for mismatches. In fact, some mismatches are actually neutral in certain ways when you consider their impact on the organism and its biology.
#
# **But how do we find the optimal solution which is by definition the longest common subsequence?** Let's not worry about that for a moment and remind ourselves of dynamic programming and recurrance relations.

# + [markdown] slideshow={"slide_type": "slide"}
# I'll put you in breakout rooms. Watch this video and then take 5 minutes to come up with a list of questions comments. Share your document in csc448 channel on Slack. I'll ask for 3-5 volunteer groups (or call on folks).
#
# <a href="https://calpoly.zoom.us/rec/share/JJuZX1D7bafCY6c07TcoHrcsOU1qA1OZ2MtEsG3qjQafOWvZKPAplPUk0i4m3kZ2.AfbmDJm0Pn2EROgL">Video from Dr. Davidson - Passcode: 7+bK*8Fv</a>

# + [markdown] slideshow={"slide_type": "slide"}
# # An Introduction to Dynamic Programming: The Change Problem
#
# **Change Problem:** Find the minimum number of coins needed to make change.
#
# Input: An integer ``money`` and an array ``coins`` of $d$ positive integers.
#
# Output: The minimum number of coins with denominations ``coins`` that changes ``money``.

# + [markdown] slideshow={"slide_type": "subslide"}
# Consider the problem to change 76 cents in a country with only three denominations: ``coins=[6,5,1]``. A minimal colection of coins totaling 76 cents must be one of the following:
# * a minimal collection of coins totaling 75 cents, plus a 1-cent coin
# * a minimal collection of coins totaling 71 cents, plust a 5-cent coin
# * a minimal collection of coins totaling 70 cents, plus a 6-cent coin
#
# This illustrates a nice recurrence relation for us:
# $$
# \mbox{MinNumCoins}(money) = min
# \left\{
#     \begin{array}{l}
#       \mbox{MinNumCoins}(money-coin_1) + 1\\
#       \mbox{...}\\
#       \mbox{MinNumCoins}(money-coin_d)+1
#     \end{array}
#   \right\}
# $$

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 3** Solve the change problem using the recurance relation given above.
#
# Input: An integer ``money`` and an array ``coins`` of $d$ positive integers.
#
# Output: The minimum number of coins with denominations ``coins`` that changes ``money``.

# + slideshow={"slide_type": "fragment"}
import numpy as np

def min_num_coins(money,coins):
    min_coins = np.Inf
    # YOUR SOLUTION HERE
    return min_coins


min_num_coins(27,[6,5,1])

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Let's time some results

# + slideshow={"slide_type": "fragment"}
# %%timeit
min_num_coins(13,[6,5,1])

# + slideshow={"slide_type": "fragment"}
# %%timeit
min_num_coins(27,[6,5,1])

# + slideshow={"slide_type": "fragment"}
# %%timeit
min_num_coins(35,[6,5,1])

# + slideshow={"slide_type": "fragment"}
# %%timeit
min_num_coins(47,[6,5,1])

# + [markdown] slideshow={"slide_type": "subslide"}
# That escalated quickly!

# + slideshow={"slide_type": "fragment"}
# %matplotlib inline
import pandas as pd
ax=pd.DataFrame({"money":[13,27,35,47],"time":[43.2e-6,5.59e-3,96.3e-3,6.24]}).plot.line(x="money",y="time")
ax.set_xlabel('Money')
ax.set_ylabel('Time');

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Is this efficient?
# Well. Of course not. It recomputes solutions to subproblems over and over again. We will solve this problem using dynamic programming. 
#
# The book takes a practical approach to dynamic programming, but I would suggest everyone read this section at least on Dynamic programming (https://en.wikipedia.org/wiki/Dynamic_programming#Computer_programming). 

# + [markdown] slideshow={"slide_type": "subslide"}
# ## Bottom-up approach to coins problem
# Let's say you knew the following solution to the coins problem for money $\le$ 12:

# + slideshow={"slide_type": "fragment"}
pd.DataFrame({"MinNumCoins(money)":[0,1,2,3,1,1,2,3,2,2,2,3,3]},index=pd.Index([0,1,2,3,4,5,6,7,8,9,10,11,12],name="money")).T


# + [markdown] slideshow={"slide_type": "subslide"}
# Could you easily compute MinNumCoins(13)? What about MinNumCoins(14)? Let's consider MinNumCoins(13). You only have three potential coins you could add at any given time: [6,5,1]. 
# * You could try to add in a 6, which would mean that MinNumCoins(13) = MinNumCoins(13-6)+1 = MinNumCoins(7)+1 = 4. 
#     * We subtract the 6 because we need to identify the optimal solution to the subproblem we need to examine. 
# * We next try 5: MinNumCoins(13) = MinNumCoins(13-5)+1 = MinNumCoins(8)+1 = 3, which is our best solution so far.
# * We next try 1: MinNumCoins(13) = MinNumCoins(13-1)+1 = MinNumCoins(12)+1 = 4, which isn't as good.
# * We don't have anything else to try, so we know our answer is 3!
# * No recursion, and a simple thing to code up.

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 4** Solve the change problem using a bottom-up dynamic programming strategy.
#
# Input: An integer ``money`` and an array ``coins`` of $d$ positive integers.
#
# Output: The minimum number of coins with denominations ``coins`` that changes ``money``.

# + slideshow={"slide_type": "fragment"}
def min_num_coins_dynamic(money,coins):
    min_coins = {0:0} # Base case, no coins needed for no money
    # YOUR SOLUTION HERE
    return min_coins[m]


min_num_coins_dynamic(27,[6,5,1])

# + [markdown] slideshow={"slide_type": "subslide"}
# Much better runtime!

# + slideshow={"slide_type": "fragment"}
# %%timeit
min_num_coins_dynamic(47,[6,5,1])

# + [markdown] slideshow={"slide_type": "slide"}
# ## Back to sequence alignment
# What is our recurence relation? Think top down. Remember our greedy actions that were available. We still only have those options available to us. For clarity, I've put our choices into a dataframe with some additional information including the part that would be added to the alignment for s1 and s2. The final column is the score that would be added to the total score.

# + slideshow={"slide_type": "fragment"}
s1 = "ACGT"
s2 = "AGCTA"
choices_df = pd.DataFrame({"remainder(s1)":[s1[1:],s1[1:],s1],
              "remainder(s2)":[s2[1:],s2,s2[1:]],
              "s1_part":[s1[0],s1[0],"-"],
              "s2_part":[s2[0],"-",s2[0]],
              "score(s1_part,s2_part)":[int(s1[0]==s2[0]),0,0]})
choices_df

# + slideshow={"slide_type": "subslide"}
s1 = "GCGT"
s2 = "AGCTA"
choices_df = pd.DataFrame({"remainder(s1)":[s1[1:],s1[1:],s1],
              "remainder(s2)":[s2[1:],s2,s2[1:]],
              "s1_part":[s1[0],s1[0],"-"],
              "s2_part":[s2[0],"-",s2[0]],
              "score(s1_part,s2_part)":[int(s1[0]==s2[0]),0,0]})
choices_df

# + [markdown] slideshow={"slide_type": "subslide"}
# Now we can write our recurrance relationship for our specific example:
#
# $$
# align(\mbox{ACGT},\mbox{AGCTA}) = max
# \left\{
#     \begin{array}{lll}
#       \mbox{align}(\mbox{CGT},\mbox{GCTA}) & + & \mbox{score}(\mbox{A},\mbox{A})\\
#       \mbox{align}(\mbox{CGT},\mbox{AGCTA}) & + & \mbox{score}(\mbox{A},\mbox{-})\\
#       \mbox{align}(\mbox{ACGT},\mbox{GCTA}) & + & \mbox{score}(\mbox{-},\mbox{A})
#     \end{array}
#   \right\}
# $$
# General case:
# $$
# align(\mbox{s1},\mbox{s2}) = max
# \left\{
#     \begin{array}{lll}
#       \mbox{align}(\mbox{s1[1:]},\mbox{s2[1:]}) & + & \mbox{score}(\mbox{s1[0]},\mbox{s2[0]})\\
#       \mbox{align}(\mbox{s1[1:]},\mbox{s2}) & + & \mbox{score}(\mbox{s1[0]},\mbox{-})\\
#       \mbox{align}(\mbox{s1},\mbox{s2[1:]}) & + & \mbox{score}(\mbox{-},\mbox{s2[0]})
#     \end{array}
#   \right\}
# $$

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 5** Solve the longest common subsequence problem using the above recurrence relation
#
# Input: Two strings
#
# Output: A longest common subsequence of these strings represented as a tuple of a (score, aligned string 1, and aligned string 2).
#
# Suggestion: If you are struggling with the recursive solution, then proceed to the dynamic programming one and come back to this.

# + slideshow={"slide_type": "fragment"}
import pandas as pd
import numpy as np
def align(s1,s2):
    # Below are the exact base cases that I want you to use
    if len(s1) == 0:
        aligned_s1 = "".join(["-" for i in range(len(s2))])
        return 0,aligned_s1,s2
    if len(s2) == 0: # no way to match
        aligned_s2 = "".join(["-" for i in range(len(s1))])
        return 0,s1,aligned_s2
    
    # You don't have to use my dataframe that helps with the choices, but ... I recommend it
    choices_df = pd.DataFrame({
        "remainder(s1)":[s1[1:],s1[1:],s1],
        "remainder(s2)":[s2[1:],s2,s2[1:]],
        "s1_part":[s1[0],s1[0],"-"],
        "s2_part":[s2[0],"-",s2[0]],
        "score(s1_part,s2_part)":[int(s1[0]==s2[0]),0,0]})
    max_score = -np.Inf
    aligned_s1 = None
    aligned_s2 = None
    for i,choice in choices_df.iterrows():
        # here is how to get these values into base Python
        rem_s1,rem_s2,s1_part,s2_part,score = choice.values
        # YOUR SOLUTION HERE
        # print(rem_s1,rem_s2)
    return max_score,aligned_s1,aligned_s2


score, aligned_s1, aligned_s2 = align("AACCT","ACACTG")
print(score)
print(aligned_s1)
print(aligned_s2)

# + [markdown] slideshow={"slide_type": "subslide"}
# ## What about a dynamic programming solution to the algorithm above?
# We can think about this the same way. We need a bottom-up approach for sequence alignment. What we need to understand the best is a good way to represent a partial solution. Consider the following matrix and what each cell might represent? What value would we like to know?

# + slideshow={"slide_type": "fragment"}
s1,s2="AACCT","ACACTG"
scores = pd.DataFrame(index=["-"]+[s1[:i+1] for i in range(len(s1))],columns=["-"]+[s2[:i+1] for i in range(len(s2))])
scores

# + [markdown] slideshow={"slide_type": "subslide"}
# Let's now fill in some values that are obvious to us.

# + slideshow={"slide_type": "fragment"}
scores.loc["-","-"] = 0
scores.loc["-","A"] = 0
scores.loc["-","AC"] = 0
# and so on... so let's do it automatically
for s2_part in scores.columns:
    scores.loc["-",s2_part] = 0
for s1_part in scores.index:
    scores.loc[s1_part,"-"] = 0
scores


# + [markdown] slideshow={"slide_type": "subslide"}
# Now how would you fill in scores.loc["A","A"]? Well... There are three options like always. These correspond to the three choices we always. We can match A to A which would obviously be good in this example. We could match the A in s1 to a gap '-' in s2. We could match A in s2 to a gap in s1. And those are our only options. If we fill out this table in a manner that moves from left to right and top to bottom, then we can figure out the max score without any issues.

# + [markdown] slideshow={"slide_type": "subslide"}
# **Exercise 6** Solve the longest common subsequence problem using dynamic programming.
#
# Input: Two strings
#
# Output: The length of the longest common subsequence of these strings.

# + slideshow={"slide_type": "subslide"}
def align_dynamic(s1,s2):
    scores = pd.DataFrame(index=["-"]+[s1[:i+1] for i in range(len(s1))],columns=["-"]+[s2[:i+1] for i in range(len(s2))])
    for s2_part in scores.columns:
        scores.loc["-",s2_part] = 0
    for s1_part in scores.index:
        scores.loc[s1_part,"-"] = 0
    
    nrows,ncols = scores.shape
    for i in range(1,nrows):
        for j in range(1,ncols):
            # What are our three options
            opt1_s1 = scores.index[i-1] # remember the rows are representative of s1
            opt1_s2 = scores.columns[j-1] # remember the columns are representative of s2
            score_opt1 = -np.Inf # FIX THIS!
            
            opt2_s1 = scores.index[i-1]
            opt2_s2 = scores.columns[j]
            score_opt2 = -np.Inf # FIT THIS!
            
            opt3_s1 = scores.index[i]
            opt3_s2 = scores.columns[j-1]
            score_opt3 = -np.Inf # FIT THIS!
            
            scores.loc[scores.index[i],scores.columns[j]] = max(score_opt1,score_opt2,score_opt3)
            
    return scores.loc[s1,s2]


score = align_dynamic("AACCT","ACACTG")
score


# + [markdown] slideshow={"slide_type": "subslide"}
# We did it! Sort of... We don't know the alignment. Only the score. We need to add this history to our algorithm.
#
# **Exercise 7** Solve the longest common subsequence problem using dynamic programming.
#
# Input: Two strings
#
# Output: A longest common subsequence of these strings represented as a tuple of a (score, aligned string 1, and aligned string 2).

# + slideshow={"slide_type": "subslide"}
def align_dynamic2(s1,s2,verbose=False):
    scores = pd.DataFrame(index=["-"]+[s1[:i+1] for i in range(len(s1))],columns=["-"]+[s2[:i+1] for i in range(len(s2))])
    aligned = pd.DataFrame(index=["-"]+[s1[:i+1] for i in range(len(s1))],columns=["-"]+[s2[:i+1] for i in range(len(s2))])
    for s2_part in scores.columns:
        scores.loc["-",s2_part] = 0
        if s2_part == "-":
            aligned.loc["-","-"] = ("","")
        else:
            aligned.loc["-",s2_part] = ("".join(["-" for i in range(len(s2_part))]),s2_part)
    for s1_part in scores.index:
        scores.loc[s1_part,"-"] = 0
        if s1_part == "-":
            aligned.loc["-","-"] = ("","")
        else:
            aligned.loc[s1_part,"-"] = (s1_part,"".join(["-" for i in range(len(s1_part))]))
    if verbose:
        display(aligned)
    
    nrows,ncols = scores.shape
    for i in range(1,nrows):
        for j in range(1,ncols):
            # What are our three options
            opt1_s1 = scores.index[i-1] # remember the rows are representative of s1
            opt1_s2 = scores.columns[j-1] # remember the columns are representative of s2
            score_opt1 = -np.Inf # FIX THIS!
            s1_aligned_opt1 = "" # FIX THIS!
            s2_aligned_opt1 = "" # FIX THIS!
            
            opt2_s1 = scores.index[i-1]
            opt2_s2 = scores.columns[j]
            score_opt2 = -np.Inf # FIT THIS!
            s1_aligned_opt2 = "" # FIX THIS!
            s2_aligned_opt2 = "" # FIX THIS!
            
            opt3_s1 = scores.index[i]
            opt3_s2 = scores.columns[j-1]
            score_opt3 = -np.Inf # FIT THIS!
            s1_aligned_opt3 = "" # FIX THIS!
            s2_aligned_opt3 = "" # FIX THIS!
            
            scores.loc[scores.index[i],scores.columns[j]] = max(score_opt1,score_opt2,score_opt3)
            if max(score_opt1,score_opt2,score_opt3) == score_opt1:
                aligned.loc[scores.index[i],scores.columns[j]] = (s1_aligned_opt1,s2_aligned_opt1)
            elif max(score_opt1,score_opt2,score_opt3) == score_opt2:
                aligned.loc[scores.index[i],scores.columns[j]] = (s1_aligned_opt2,s2_aligned_opt2)
            else:
                aligned.loc[scores.index[i],scores.columns[j]] = (s1_aligned_opt3,s2_aligned_opt3)
    if verbose:
        display(scores)
        display(aligned)
    return scores.loc[s1,s2],aligned.loc[s1,s2][0],aligned.loc[s1,s2][1]


score,s1_aligned,s2_aligned = align_dynamic2("AACCT","ACACTG")
print(score)
print(s1_aligned)
print(s2_aligned)

# + [markdown] slideshow={"slide_type": "subslide"}
# <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/0/02/BLOSUM62.png/400px-BLOSUM62.png" width=800>

# + slideshow={"slide_type": "skip"}
# Don't forget to push!
# -


