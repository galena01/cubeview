// Initialize 3Dmol.js viewer
const viewer = $3Dmol.createViewer(
  document.getElementById("viewer"),
  { backgroundColor: "white" }
);

async function loadMolecule() {
  try {
    const response = await fetch("./mol.xyz");
    if (!response.ok) throw new Error("Failed to load molecule structure");
    const xyzData = await response.text();

    const model = viewer.addModel(xyzData, "xyz");
    viewer.addStyle({}, { stick: { radius: 0.05 } });

    const atoms = model.selectedAtoms({});
    let index = 0;
    for (let atom of atoms) {
      const labelText = index + " " + atom.elem;
      index += 1;

      if (show_atom_index) {
        viewer.addLabel(labelText, {
          position: { x: atom.x, y: atom.y, z: atom.z },
          fontColor: 'black',
          fontSize: 16,
          inFront: true,
          alignment: 'center',
          showBackground: false
        });
      }
    }

    viewer.zoomTo();
    viewer.render();
    return true;
  } catch (error) {
    console.error("Error loading molecule structure:", error);
    return false;
  }
}

async function loadCube(spin, index) {
  try {
    viewer.clear();
    document.getElementById("error-message").style.display = "none";

    const molOk = await loadMolecule();
    if (!molOk) throw new Error("Molecule structure loading failed");

    let fileName = "";
    if (spin && spin !== "") {
      fileName = `${spin}_${index}.cube.gz`;
    } else {
      fileName = `${index}.cube.gz`;
    }

    if (spin && spin !== "") {
      document.getElementById("info-box").textContent =
        `Spin: ${spin} | Index: ${index}`;
    } else {
      document.getElementById("info-box").textContent =
        `Index: ${index}`;
    }

    const mask = document.getElementById("loadingMask");
    mask.style.visibility = "visible";

    const response = await fetch(`./cubes/${fileName}`);
    if (!response.ok) {
      if (response.status === 404) {
        document.getElementById("error-message").textContent =
          `Cube file not rendered: ${fileName}`;
        document.getElementById("error-message").style.display = "block";
      } else {
        throw new Error(`File loading failed: ${response.status}`);
      }
    }
    const arrayBuffer = await response.arrayBuffer();
    const cubeData = pako.inflate(new Uint8Array(arrayBuffer), { to: 'string' });

    viewer.addVolumetricData(cubeData, "cube", {
      isoval: isovalue,
      color: "red",
      opacity: 0.75,
    });
    viewer.addVolumetricData(cubeData, "cube", {
      isoval: -isovalue,
      color: "blue",
      opacity: 0.75,
    });

    viewer.zoomTo();
    viewer.render();
  } catch (error) {
    console.error("Error loading cube file:", error);
  } finally {
    document.getElementById("loadingMask").style.visibility = "hidden";
  }
}

async function loadOrbitalsJSON() {
  try {
    const response = await fetch("orbitals.json");
    if (!response.ok) throw new Error("Failed to load orbitals.json");
    const data = await response.json();

    data.sort((a, b) => {
      const spinA = (a.spin || "").toLowerCase();
      const spinB = (b.spin || "").toLowerCase();
      if (spinA < spinB) return -1;
      if (spinA > spinB) return 1;
      return parseInt(a.index, 10) - parseInt(b.index, 10);
    });

    let tableHTML = `
    <div class="table-container">
      <table>
        <thead>
          <tr>
            <th>MO #</th>
            <th>Spin</th>
            <th>Irrep</th>
            <th>Occ</th>
            <th>Energy</th>
            ${data.length > 0 && 'ao_components' in data[0] ? '<th>AO Components</th>' : ''}
          </tr>
        </thead>
        <tbody>
    `;
    data.forEach((orbital) => {
      const spinText = orbital.spin ? orbital.spin : "a/b";
      tableHTML += `
      <tr data-orbital='${JSON.stringify(orbital)}'>
        <td>${orbital.index}</td>
        <td>${spinText}</td>
        <td>${orbital.irrep}</td>
        <td>${parseFloat(orbital.occ).toFixed(4)}</td>
        <td>${parseFloat(orbital.energy).toFixed(4)}</td>
        ${'ao_components' in orbital ? `<td>${orbital.ao_components}</td>` : ''}
      </tr>
      `;
    });
    tableHTML += `
        </tbody>
      </table>
    </div>
  `;

    document.getElementById("sidebar").innerHTML = tableHTML;

    const rows = document.querySelectorAll("#sidebar table tbody tr");
    rows.forEach((row) => {
      row.addEventListener("click", () => {
        rows.forEach((r) => r.classList.remove("selected"));
        row.classList.add("selected");

        const orbitalData = JSON.parse(row.getAttribute('data-orbital'));
        const idx = orbitalData.index;
        const spin = orbitalData.spin || '';
        
        const spinType = spin.toLowerCase() === "b" ? "beta" :
          spin.toLowerCase() === "a" ? "alpha" :
          spin.toLowerCase() === "a/b" ? "" : spin; 
        loadCube(spinType, idx);

      });
    });
  } catch (err) {
    document.getElementById("sidebar").innerHTML =
      `<p style="color:red;">Error: ${err.message}</p>`;
    console.error(err);
  }
}

let isovalue = 0.05;
let show_atom_index = false; 

const isoSlider = document.getElementById('isovalue');
const isoLabel = document.getElementById('isovalue-label');
isoSlider.addEventListener('input', () => {
  isovalue = parseFloat(isoSlider.value);
  isoLabel.textContent = isovalue.toFixed(2);
  
  const selectedRow = document.querySelector('#sidebar table tbody tr.selected');
  if (selectedRow) {
    const orbitalData = JSON.parse(selectedRow.getAttribute('data-orbital'));
    const idx = orbitalData.index;
    const spin = orbitalData.spin || '';
    const spinType = spin.toLowerCase() === "b" ? "beta" :
      spin.toLowerCase() === "a" ? "alpha" :
      spin.toLowerCase() === "a/b" ? "" : spin; 
    loadCube(spinType, idx);
  }
});

const showAtomIndexCheckbox = document.getElementById('show-atom-index');
showAtomIndexCheckbox.addEventListener('change', function () {
  const checked = showAtomIndexCheckbox.checked;
  if (checked) {
    show_atom_index = true;
  } else {
    show_atom_index = false;
  }

  const selectedRow = document.querySelector('#sidebar table tbody tr.selected');
  if (selectedRow) {
    const orbitalData = JSON.parse(selectedRow.getAttribute('data-orbital'));
    const idx = orbitalData.index;
    const spin = orbitalData.spin || '';
    const spinType = spin.toLowerCase() === "b" ? "beta" :
      spin.toLowerCase() === "a" ? "alpha" :
      spin.toLowerCase() === "a/b" ? "" : spin; 
    loadCube(spinType, idx);
  }
});

// Page load immediately calls
document.addEventListener("DOMContentLoaded", () => {
  loadOrbitalsJSON();
});