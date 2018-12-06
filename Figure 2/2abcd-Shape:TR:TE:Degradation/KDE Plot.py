
#############
#  mouse
#############

### GINI - TE

scatter_plot(mouse_5UTR_gini['cy'], mouse_TE, color="#8172B2", title="", xlabel="cy_5UTR_gini", ylabel="TE", xlim=[0.4, 1.0], ylim=[-2, 2], Parser=mm10_parseTrans, kind='kde')
plt.savefig("figs/mouse_TE_kde.pdf")
plt.show()

### GINI - TS

scatter_plot(mouse_5UTR_gini['ch'], mouse_TS, color="#55A868", title="", xlabel="ch_5UTR_gini", ylabel="TS", xlim=[0.4, 1.0], ylim=[0, 10], Parser=mm10_parseTrans, kind='kde')
plt.savefig("figs/mouse_TS_kde.pdf")
plt.show()


#############
#  human
#############

### GINI - degaradation

scatter_plot(human_gini['np'], human_hl, color="#C44E52", title="", xlabel="np_gini", ylabel="log2(half-life)", xlim=[0.4, 1.0], ylim=[4, 10], Parser=hg38_parseTrans, ylog=True, kind='kde')
plt.savefig("figs/human_HL_np_kde.pdf")
plt.show()

scatter_plot(human_gini['cy'], human_hl, color="#8172B2", title="", xlabel="cy_gini", ylabel="log2(half-life)", xlim=[0.4, 1.0], ylim=[4, 10], Parser=hg38_parseTrans, ylog=True, kind='kde')
plt.savefig("figs/human_HL_cy_kde.pdf")
plt.show()

