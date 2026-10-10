function savefigs(fig, fig_folder, name)

saveas(fig, fig_folder + "\" + name, "pdf");
saveas(gcf, fig_folder + "\" + name, "png");
saveas(fig, fig_folder + "\" + name, "eps");

end