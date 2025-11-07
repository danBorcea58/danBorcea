


def toroydConnections(STL, GEO, CH, np):
    steps = 3
    STL["interconnection"] = {
        "v": np.empty((0, 3), dtype=float),
        "f": np.empty((0, 3), dtype=int)
    }

    for i in range(int(CH["number"] - 1)):
        start_prev = STL["toroydConnection"][i]["end"]          # (N,3)
        start_curr = STL["toroydConnection"][i]["start"]        # (N,3)
        centreStart = STL["toroydConnection"][i]["start_centre3D"]  # (N,3)
        N = start_prev.shape[0]

        # direzione verso il centro di curvatura
        dir_to_centre = centreStart - start_curr
        norm_dir = np.linalg.norm(dir_to_centre, axis=1)[:, None]
        norm_dir[norm_dir == 0] = 1
        dir_to_centre /= norm_dir  # normalizza

        # distanza media e passo iniziale
        diff_init = start_prev - start_curr
        mean_dist = np.mean(np.linalg.norm(diff_init, axis=1))
        centre_step = mean_dist / (2 * steps)

        # Genera i layer intermedi adattivi
        all_layers = [start_curr.copy()]
        layer = start_curr.copy()

        for k in range(1, 2 * steps + 1):
            if k == 1:
                # spostamento verso il centro di curvatura
                layer = layer + dir_to_centre * centre_step
            else:
                # aggiorna direzione e distanza residua verso il target
                remaining = start_prev - layer
                step_vec = remaining / (2 * steps - k + 1)  # passo adattivo
                layer = layer + step_vec

            all_layers.append(layer.copy())

        # aggiorna i vertici globali
        new_vertices = np.vstack(all_layers)
        base_idx = len(STL["interconnection"]["v"])
        STL["interconnection"]["v"] = (
            np.vstack([STL["interconnection"]["v"], new_vertices])
            if STL["interconnection"]["v"].size
            else new_vertices
        )

        # Costruisci le facce tra anelli consecutivi
        faces = []
        for k in range(1, len(all_layers)):
            curr_idx = np.arange(base_idx + k*N, base_idx + (k+1)*N)
            prev_idx = np.arange(base_idx + (k-1)*N, base_idx + k*N)
            for j in range(N - 1):
                faces.append([curr_idx[j], prev_idx[j], prev_idx[j+1]])
                faces.append([curr_idx[j], prev_idx[j+1], curr_idx[j+1]])


        STL["interconnection"]["f"] = np.vstack([STL["interconnection"]["f"], np.array(faces, dtype=int)])

    return STL