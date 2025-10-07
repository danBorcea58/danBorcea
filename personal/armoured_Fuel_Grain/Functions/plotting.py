



def plotFunc(SEC, GEO, PARAM, np, plt):
    SEC["plotNum"]          = {}
    SEC["plotX"]            = {}
    SEC["plotY"]            = {}

    X_all = []
    Y_all = []

    for i in range(SEC["number"]):
        SEC["plotNum"][i] = round(SEC["length"][i] / (GEO["radius"]-GEO["port"])*500)
        SEC["plotX"][i] = np.linspace(SEC["port"][i], SEC["radius"][i], SEC["plotNum"][i])
        SEC["plotY"][i] = SEC["waveLength"][i] - (SEC["plotX"][i] - SEC["port"][i]) * SEC["rate"][i]

        # concatena i punti
        X_all.append(SEC["plotX"][i])
        Y_all.append(SEC["plotY"][i])

    # Unisci in un unico array
    X_all = np.concatenate(X_all)
    Y_all = np.concatenate(Y_all)

    # Regressione polinomiale di terzo grado
    coeffs = np.polyfit(X_all, Y_all, 10)  
    PARAM["coeffs"] = coeffs 
    poly_fit = np.poly1d(coeffs)

    # Valutazione della curva
    X_fit = np.linspace(min(X_all), max(X_all), 500)
    Y_fit = poly_fit(X_fit)

    plt.figure(figsize=(8,5))
    plt.plot(X_fit, 
            Y_fit, 
            '-',
            color = "salmon", 
            linewidth = 3, 
            label="Polynomial fit (10th order)")
    plt.plot(X_all, 
            Y_all, 
            ':', 
            color = "blue",
            markersize = 1, 
            label="Calibration trace")
    plt.xlabel("Radial position [mm]")
    plt.ylabel("Wavelength [mm]")
    plt.title("Gyroid wavelength variation along radius")
    plt.grid(True)
    plt.legend()

    plt.savefig("./Output/Chart/gyroid_wavelength_fit.png", dpi=300, bbox_inches="tight")
    plt.close()
    return PARAM