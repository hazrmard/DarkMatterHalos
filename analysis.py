from halos import *
from halos.helpers import *
#from scipy import stats


def eval_vs_cleave(H,order=1):
    """
    takes a group of halos and checks how well radii ratios based on eigenvalues
    fit with ratios based on maximum dimensions. Assumes halos are NOT self-similar
    i.e. transform=True.
    eigenvalue ratio = ratio of median and max evals along principal axis
    dimension ratio = ratio of x dimension of inner halo to all halo particles
    :H a HalfMassRadius instance
    :order order of fitting
    """
    H.higher_order_fit(order=order)
    e_ratio = [h.half_mass_radius/max(h.radii) for h in H.halos]
    d_ratio = [h.inner_R.x/h.cleave().x for h in H.halos]
    corr_coeff = np.corrcoef(e_ratio, d_ratio)
    plt.scatter(e_ratio, d_ratio, alpha=0.3)
    plt.xlabel('Eigenvalue ratios')
    plt.ylabel('Dimension ratios')
    plt.suptitle('Correlation Coefficient: ' + '{:.3f}'.format(corr_coeff[0,1]))
    plt.title('Eigenvalue vs. Cleave relationship')
    plt.show()
    plt.close()


def radius_distribution(H, mode='eval', bins=10):
    """Calculate the distribution of inner to outer radius ratios for halos.
    Ratios can be computed based on eigenvalues or absolute dimensions.
    :H a HalfMassRadius object containing halos in H.halos
    :mode 'eval' for eigenvalue ratios, 'cleave' for absolute dimension ratios
    :bins number of uniformly distributed bins
    """
    if mode=='eval':
        ratios = [h.half_mass_radius/max(h.radii) for h in H.halos]
    elif mode=='cleave':
        ratios = [h.inner_R.x/h.cleave().x for h in H.halos]
    else:
        return -1
    counts, bin_edges = np.histogram(ratios, bins)
    rinds = np.digitize(ratios, bin_edges)
    rinds = [r-1 for r in rinds] #identify bins by lower bound
    rsum = [0]*(len(counts))
    rssum = [0]*(len(counts))
    rmean = [0]*(len(counts))
    rstd_dev = [0]*(len(counts))
    m = max(rinds)
    rinds = [x-1 if x==m else x for x in rinds] # correction to make last bin inclusive of right edge
    for i in range(len(rinds)):
        rsum[rinds[i]] += ratios[i]
        rssum[rinds[i]] += ratios[i]**2
    for i in range(len(counts)):
        if counts[i]!=0:
            rmean[i] = rsum[i] / counts[i]
            rstd_dev[i] = np.sqrt(rssum[i]/counts[i] - rmean[i]**2 )
    bins_mean = [0.5*(bin_edges[i]+bin_edges[i+1]) for i in range(len(bin_edges)-1)]
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    #ax1.scatter(rmean, counts)
    arg_max = np.argmax(counts)
    ax1.set_title('Half Mass Radii - Max:(' + str(counts[arg_max]) + ',' + '{:.3f}'.format(rmean[arg_max]) + '); Total:' + str(len(H.halos)))
    ax1.errorbar(rmean, counts, xerr=rstd_dev, fmt='o')
    ax1.set_xlabel('Inner to Outer halo radius')
    ax1.set_ylabel('Frequency')
    fig.show()
    fig.close()


def particles_vs_radii(H, mode='eval'):
    """Plot the relationship between halo size and its inner:outer halo radius.
    :H HalfMassRadius object with halo objects in H.halos
    :mode either to use eigenvalues or absolute dimensions for calculating ratios.
    """
    particles = [len(h.particles) for h in H.halos]
    if mode=='eval':
        ratios = [h.half_mass_radius/max(h.radii) for h in H.halos]
    elif mode=='cleave':
        ratios = [h.inner_R.x/h.cleave().x for h in H.halos]
    else:
        return -1
    corr_coeff = np.corrcoef(ratios, particles)
    plt.scatter(ratios, particles, alpha=0.3)
    plt.xlabel('Radii (Inner to Outer)')
    plt.ylabel('# of halo particles')
    plt.suptitle('Correlation Coefficient: ' + '{:.3f}'.format(corr_coeff[0,1]))
    plt.title('Halo Size vs. Radius ratios')
    plt.show()
    plt.close()
