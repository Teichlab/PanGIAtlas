import networkx as nx
import numpy as np
import plotnine as p9
import matplotlib.pyplot as plt
from adjustText import adjust_text


def plot_network(cell_pairs,
                 score_key,
                 top_n=20, 
                 figsize=(7, 7),
                 network_layout='kamada_kawai', 
                 edge_alpha=0.5, 
                 edge_arrow_size=10,
                 edge_width_factor=1,
                 edge_color='Blues',
                 node_color='skyblue', 
                 node_size=350, 
                 node_alpha=0.9, 
                 node_label_alpha=0.8,
                 node_label_size=9, 
                 node_label_offset=[0.0, 0.15], 
                 title_size=14, 
                 centralize_nodes=None, 
                 filename=None,
                 show_colorbar=False,
                 scale_weights=True,
                 **kwargs,
                 ):
    
    fig, ax = plt.subplots(figsize=figsize)
    
    cp = cell_pairs.copy()
    if scale_weights:
        cp[score_key] = (cp[score_key] - cp[score_key].min()) / (cp[score_key].max() - cp[score_key].min()) + 1
    cp = cp.sort_values(score_key, ascending=False).head(top_n)
    
    cp['source'] = cp['source'].astype('category')
    cp['target'] = cp['target'].astype('category')

    G = nx.from_pandas_edgelist(cp, source='source', target='target', edge_attr=score_key, create_using=nx.DiGraph())

    # Determine layout
    if network_layout == 'spring':
        pos = nx.spring_layout(G, **kwargs)
    elif network_layout == 'circular':
        pos = nx.circular_layout(G, **kwargs)
    elif network_layout == 'spectral':
        pos = nx.spectral_layout(G, **kwargs)
    elif network_layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G, **kwargs)
    else:
        raise ValueError("network_layout should be either 'spring' or 'circular'")
    
    # Adjust positions if centralization is needed
    if centralize_nodes in ['source', 'target']:
        selected_nodes = cp[centralize_nodes].unique()
        for node in selected_nodes:
            pos[node] = np.mean([pos[n] for n in selected_nodes], axis=0) + np.random.rand(2) * 0.5

    # Scale the edge weights
    edge_weights = np.array([G[u][v][score_key] for u, v in G.edges()]) * edge_width_factor

    # Create a colormap
    cmap = plt.cm.get_cmap(edge_color)
    
    # Draw network with directed edges
    nx.draw_networkx_edges(G, pos, alpha=edge_alpha, arrows=True, arrowsize=edge_arrow_size, width=edge_weights,
                           edge_color=edge_weights, edge_cmap=cmap, ax=ax, connectionstyle="arc3,rad=-0.25",
                           edge_vmax=cp[score_key].max(), edge_vmin=cp[score_key].min())
    nx.draw_networkx_nodes(G, pos, node_color=node_color, node_size=node_size, alpha=node_alpha, ax=ax)

    label_options = {"ec": "k", "fc": "white", "alpha": node_label_alpha}
    nx.draw_networkx_labels(G, {k: v + np.array(node_label_offset) for k, v in pos.items()},
                            font_size=node_label_size, bbox=label_options, ax=ax)

    # Adjust axis limits
    ax.set_frame_on(False)
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    coeff = 1.4
    ax.set_xlim((xlim[0] * coeff, xlim[1] * coeff))
    ax.set_ylim((ylim[0] * coeff, ylim[1] * coeff))
    ax.set_title(score_key, fontsize=title_size, fontweight='bold')

    # show legend or not
    if show_colorbar:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=edge_weights.min(), vmax=edge_weights.max()))
        plt.colorbar(sm, ax=ax, label='Edge Weight', shrink=0.8)

    plt.tight_layout()
    
    if filename is not None:
        plt.savefig(filename, dpi=200, bbox_inches='tight')
    
    return fig, ax



def plot_lr_pairs(H, fct, label_fun, method='average', filename=None):
    lr_pairs = H[fct].reset_index()
    lr_pairs['rank'] = lr_pairs[fct].rank(ascending=False, method=method)
    lr_pairs['name'] = lr_pairs.apply(label_fun, axis=1)

    # Create scatter plot
    fig, ax = plt.subplots(figsize=(6, 6))  # Increased figure size for better readability
    ax.scatter(lr_pairs['rank'], lr_pairs[fct], s=30, c='black')

    # Add text labels only for non-None labels
    texts = []
    for i, row in lr_pairs.iterrows():
        if row['name']:  # only add a text label if the name is not None
            texts.append(ax.text(row['rank'], row[fct], row['name'], fontsize=14))  # Increased font size

    # Adjust text to avoid overlap with improved parameters
    adjust_text(texts, expand_points=(2, 1), arrowprops=dict(arrowstyle='->', color='darkred'), force_text=(0.5, 1))

    # Set labels and title with larger font size
    ax.set_xlabel('Rank', fontsize=15)
    ax.set_ylabel(f'{fct} loadings', fontsize=15)
    plt.xticks(fontsize=14, rotation=90)  # Rotate and increase font size of x-axis labels
    plt.yticks(fontsize=14)  # Increase font size of y-axis labels
    plt.tight_layout()  # Adjust the padding between and around subplots
    
    
    if filename:
        plt.savefig(filename, bbox_inches='tight')

    plt.show()