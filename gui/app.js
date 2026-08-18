const views = document.querySelectorAll(".view");
const navItems = document.querySelectorAll(".nav-item");

function switchView(name) {
  views.forEach(view => view.classList.toggle("active", view.id === `${name}-view`));
  navItems.forEach(item => item.classList.toggle("active", item.dataset.view === name));
  window.scrollTo({ top: 0, behavior: "smooth" });
}

navItems.forEach(item => item.addEventListener("click", () => switchView(item.dataset.view)));
document.querySelectorAll("[data-target]").forEach(button => {
  button.addEventListener("click", () => switchView(button.dataset.target));
});

function openModal(modal) {
  modal.classList.add("open");
  modal.setAttribute("aria-hidden", "false");
}

function closeModal(modal) {
  modal.classList.remove("open");
  modal.setAttribute("aria-hidden", "true");
}

const analysisModal = document.getElementById("analysisModal");
const systemModal = document.getElementById("systemModal");
document.getElementById("newAnalysisButton").addEventListener("click", () => openModal(analysisModal));
document.querySelectorAll(".module-card").forEach(card => {
  card.addEventListener("click", () => openModal(analysisModal));
});
document.querySelectorAll("#addSystemButton, #addSystemCard").forEach(button => {
  button.addEventListener("click", () => openModal(systemModal));
});
document.querySelectorAll(".modal-close").forEach(button => {
  button.addEventListener("click", () => closeModal(button.closest(".modal-backdrop")));
});
document.querySelectorAll(".modal-backdrop").forEach(backdrop => {
  backdrop.addEventListener("click", event => {
    if (event.target === backdrop) closeModal(backdrop);
  });
});
document.querySelectorAll(".modal-options button").forEach(button => {
  button.addEventListener("click", () => {
    closeModal(analysisModal);
    switchView("analyses");
    showToast("Workflow selected", `${button.dataset.module} is ready to configure.`);
  });
});
const drawer = document.getElementById("windowDrawer");
const drawerShade = document.getElementById("drawerShade");
function toggleDrawer(open) {
  drawer.classList.toggle("open", open);
  drawerShade.classList.toggle("open", open);
}
document.getElementById("editWindowButton").addEventListener("click", () => toggleDrawer(true));
drawer.querySelector(".modal-close").addEventListener("click", () => toggleDrawer(false));
drawerShade.addEventListener("click", () => toggleDrawer(false));
drawer.querySelector(".drawer-save").addEventListener("click", () => {
  toggleDrawer(false);
  showToast("Sampling policy updated", `Final ${windowRange.value} ns will be used for every replica.`);
});
const windowRange = document.getElementById("windowRange");
windowRange.addEventListener("input", () => {
  document.getElementById("windowValue").textContent = `${windowRange.value} ns`;
});

const toast = document.getElementById("toast");
let toastTimer;
function showToast(title, detail) {
  toast.querySelector("strong").textContent = title;
  toast.querySelector("small").textContent = detail;
  toast.classList.add("show");
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => toast.classList.remove("show"), 3300);
}
document.getElementById("validateButton").addEventListener("click", async event => {
  const button = event.currentTarget;
  const original = button.innerHTML;
  button.innerHTML = `<span class="queue-state running"></span> Auditing replicas`;
  try {
    await loadSystems();
    showToast("Inputs validated", `${localSystems.length} local systems were re-audited.`);
  } catch (error) {
    showToast("Local engine unavailable", error.message);
  } finally {
    button.innerHTML = original;
  }
});

document.querySelectorAll(".check-card").forEach(card => {
  card.addEventListener("click", () => {
    const input = card.querySelector("input");
    setTimeout(() => card.classList.toggle("selected", input.checked), 0);
  });
});
document.querySelectorAll(".analysis-option").forEach(option => {
  option.addEventListener("click", () => {
    const input = option.querySelector("input");
    setTimeout(() => {
      option.classList.toggle("selected", input.checked);
      const count = document.querySelectorAll(".analysis-option input:checked").length;
      document.getElementById("moduleCount").textContent = count;
      document.getElementById("runWorkflowButton").childNodes[0].textContent = `Run ${count} modules `;
    }, 0);
  });
});

function startWorkflow() {
  switchView("jobs");
  showToast("Analysis queued", "36 replica-level tasks were created.");
}
document.querySelectorAll("#runWorkflowButton, #runWorkflowSideButton").forEach(button => {
  button.addEventListener("click", startWorkflow);
});

document.addEventListener("keydown", event => {
  if (event.key === "Escape") {
    document.querySelectorAll(".modal-backdrop.open").forEach(closeModal);
    toggleDrawer(false);
  }
});

let localSystems = [];
let engineOnline = false;
let jobPoller;

async function api(path, options = {}) {
  const response = await fetch(path, {
    headers: { "Content-Type": "application/json", ...(options.headers || {}) },
    ...options
  });
  const contentType = response.headers.get("content-type") || "";
  if (!contentType.includes("application/json")) {
    throw new Error("Open the GUI through server.py to use the local analysis engine.");
  }
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || `Request failed (${response.status})`);
  return data;
}

function escapeHtml(value) {
  return String(value ?? "").replace(/[&<>"']/g, character => ({
    "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"
  })[character]);
}

async function checkEngine() {
  const badge = document.getElementById("localEngineBadge");
  const statusText = document.getElementById("engineStatus");
  try {
    const status = await api("/api/status");
    engineOnline = true;
    statusText.textContent = status.vmd_available ? "Local runner + VMD ready" : "Server ready · VMD missing";
    badge.textContent = status.vmd_available ? "Engine + VMD ready" : "Engine ready · VMD missing";
    badge.className = "engine-badge ready";
    return status;
  } catch (error) {
    engineOnline = false;
    statusText.textContent = "Frontend demonstration";
    badge.textContent = "Local engine offline";
    badge.className = "engine-badge offline";
    return null;
  }
}

function renderSystems() {
  const list = document.getElementById("localSystemList");
  const select = document.getElementById("rmsdSystem");
  if (!localSystems.length) {
    list.innerHTML = `<div class="empty-state"><span>＋</span><strong>No local systems configured</strong><small>Use “Add system” to validate topology and trajectories.</small></div>`;
    select.innerHTML = `<option value="">Import a system first</option>`;
    renderReplicaChoices();
    return;
  }
  list.innerHTML = localSystems.map(system => {
    const totalSegments = system.replicas.reduce((sum, replica) => sum + replica.segment_count, 0);
    const warnings = system.replicas.reduce((sum, replica) => sum + replica.warnings, 0);
    const replicas = system.replicas.map(replica =>
      `<span>${escapeHtml(replica.name)} · ${replica.segment_count} DCD${replica.warnings ? " · ⚠" : ""}</span>`
    ).join("");
    return `<div class="local-system-row">
      <div><strong>${escapeHtml(system.name)}</strong><small>${escapeHtml(system.description || "Local MD system")}</small></div>
      <div><strong>${system.replicas.length} replicas · ${totalSegments} segments</strong><small>${escapeHtml(system.psf)}</small></div>
      <div class="replica-summary">${replicas}</div>
      <b class="health ${warnings ? "warning" : "ready"}">${warnings ? `${warnings} warning` : "Validated"}</b>
    </div>`;
  }).join("");
  select.innerHTML = `<option value="">Choose a local system</option>` + localSystems.map(system =>
    `<option value="${escapeHtml(system.id)}">${escapeHtml(system.name)} · ${system.replicas.length} replicas</option>`
  ).join("");
  renderReplicaChoices();
}

async function loadSystems() {
  if (!engineOnline) await checkEngine();
  if (!engineOnline) {
    renderSystems();
    return;
  }
  const data = await api("/api/systems");
  localSystems = data.systems || [];
  renderSystems();
}

function importPayload() {
  const replicas = [...document.querySelectorAll(".replica-input-block")].map(block => ({
    name: block.querySelector('[data-field="name"]').value.trim(),
    trajectory_dir: block.querySelector('[data-field="trajectory_dir"]').value.trim(),
    pattern: block.querySelector('[data-field="pattern"]').value.trim(),
    frame_interval_ns: Number(block.querySelector('[data-field="frame_interval_ns"]').value),
    start_time_ns: Number(block.querySelector('[data-field="start_time_ns"]').value)
  }));
  return {
    name: document.getElementById("importSystemName").value.trim(),
    psf: document.getElementById("importPsf").value.trim(),
    pdb: document.getElementById("importPdb").value.trim(),
    replicas
  };
}

function renumberReplicas() {
  document.querySelectorAll(".replica-input-block").forEach((block, index) => {
    block.querySelector(".form-divider span").textContent = `Replica ${index + 1}`;
    const name = block.querySelector('[data-field="name"]');
    if (!name.value || /^Run \d+$/.test(name.value)) name.value = `Run ${index + 1}`;
    const remove = block.querySelector(".remove-replica");
    if (remove) remove.hidden = index === 0 && document.querySelectorAll(".replica-input-block").length === 1;
  });
}

document.getElementById("addReplicaButton").addEventListener("click", () => {
  const container = document.getElementById("replicaBlocks");
  const source = container.querySelector(".replica-input-block");
  const clone = source.cloneNode(true);
  clone.querySelectorAll("input").forEach(input => {
    if (input.dataset.field === "pattern") input.value = "md-*.dcd";
    else if (input.dataset.field === "frame_interval_ns") input.value = "0.02";
    else if (input.dataset.field === "start_time_ns") input.value = "0";
    else input.value = "";
  });
  const remove = document.createElement("button");
  remove.type = "button";
  remove.className = "remove-replica";
  remove.textContent = "×";
  remove.addEventListener("click", () => {
    clone.remove();
    renumberReplicas();
  });
  clone.appendChild(remove);
  container.appendChild(clone);
  renumberReplicas();
});

function importFeedback(message, type) {
  const feedback = document.getElementById("importFeedback");
  feedback.textContent = message;
  feedback.className = `import-feedback show ${type}`;
}

document.getElementById("validateSystemButton").addEventListener("click", async event => {
  const button = event.currentTarget;
  button.disabled = true;
  button.textContent = "Reading trajectories…";
  try {
    const data = await api("/api/systems/validate", {
      method: "POST", body: JSON.stringify(importPayload())
    });
    const replicas = data.system.replicas;
    const segments = replicas.reduce((sum, replica) => sum + replica.segment_count, 0);
    const frames = replicas.reduce((sum, replica) => sum + (replica.total_frames || 0), 0);
    const warnings = replicas.reduce((sum, replica) => sum + replica.warnings, 0);
    importFeedback(
      `Validated ${replicas.length} replica(s) and ${segments} DCD files` +
      `${frames ? ` (${frames.toLocaleString()} frames)` : ""}` +
      `${warnings ? `; ${warnings} unusually small file(s)` : "; no size anomalies"}.`,
      "success"
    );
  } catch (error) {
    importFeedback(error.message, "error");
  } finally {
    button.disabled = false;
    button.textContent = "Validate paths";
  }
});

document.getElementById("saveSystemButton").addEventListener("click", async event => {
  const button = event.currentTarget;
  button.disabled = true;
  button.textContent = "Validating and saving…";
  try {
    const data = await api("/api/systems", {
      method: "POST", body: JSON.stringify(importPayload())
    });
    closeModal(systemModal);
    await loadSystems();
    showToast("System saved", `${data.system.name} is ready for local analysis.`);
  } catch (error) {
    importFeedback(error.message, "error");
  } finally {
    button.disabled = false;
    button.innerHTML = `Save system <span>→</span>`;
  }
});

function renderReplicaChoices() {
  const systemId = document.getElementById("rmsdSystem").value;
  const system = localSystems.find(item => item.id === systemId);
  const container = document.getElementById("rmsdReplicaChoices");
  if (!system) {
    container.innerHTML = `<small>Select a local system to choose replicas.</small>`;
    return;
  }
  container.innerHTML = system.replicas.map(replica =>
    `<label><input type="checkbox" value="${escapeHtml(replica.id)}" checked> ${escapeHtml(replica.name)} · ${replica.segment_count} segments</label>`
  ).join("");
}
document.getElementById("rmsdSystem").addEventListener("change", renderReplicaChoices);

document.getElementById("runLocalRmsdButton").addEventListener("click", async event => {
  const button = event.currentTarget;
  const systemId = document.getElementById("rmsdSystem").value;
  const replicaIds = [...document.querySelectorAll("#rmsdReplicaChoices input:checked")].map(input => input.value);
  if (!systemId) return showToast("Choose a local system", "Import or select a validated system first.");
  button.disabled = true;
  button.textContent = "Creating RMSD job…";
  try {
    const selection = document.getElementById("rmsdSelection").value.trim();
    const data = await api("/api/rmsd/run", {
      method: "POST",
      body: JSON.stringify({
        system_id: systemId,
        replica_ids: replicaIds,
        config: {
          stride: Number(document.getElementById("rmsdStride").value) || 1,
          alignment_selection: document.getElementById("rmsdAlignment").value.trim(),
          selections: [{ name: "rmsd_A", selection }]
        }
      })
    });
    switchView("jobs");
    await loadJobs();
    showToast("RMSD job started", `${data.job.system_name} is running through VMD.`);
  } catch (error) {
    showToast("RMSD job not started", error.message);
  } finally {
    button.disabled = false;
    button.innerHTML = `Run RMSD locally <span>→</span>`;
  }
});

function jobDuration(job) {
  const end = job.completed_at ? new Date(job.completed_at) : new Date();
  const start = new Date(job.created_at);
  const seconds = Math.max(0, Math.round((end - start) / 1000));
  return `${String(Math.floor(seconds / 60)).padStart(2, "0")}:${String(seconds % 60).padStart(2, "0")}`;
}

function renderJobs(jobs) {
  const list = document.getElementById("localJobsList");
  if (!jobs.length) {
    list.innerHTML = `<div class="job-empty">No local jobs yet. Import a system and run RMSD.</div>`;
    return;
  }
  list.innerHTML = jobs.map(job => {
    const state = job.status === "completed" ? "ready" : job.status === "failed" ? "warning" : "running";
    return `<div class="job-line">
      <strong>${escapeHtml(job.type.toUpperCase())}</strong>
      <span>${escapeHtml(job.system_name)} · ${job.replica_ids.length} replica(s)</span>
      <div class="job-progress" title="${escapeHtml(job.message)}"><i style="width:${job.progress || 0}%"></i></div>
      <span>${jobDuration(job)}</span>
      <b class="health ${state}">${escapeHtml(job.status)}</b>
    </div>`;
  }).join("");
}

async function loadJobs() {
  if (!engineOnline) return;
  const data = await api("/api/jobs");
  renderJobs(data.jobs || []);
  const active = (data.jobs || []).some(job => ["queued", "running"].includes(job.status));
  clearTimeout(jobPoller);
  if (active) jobPoller = setTimeout(loadJobs, 2000);
}

async function initializeLocalEngine() {
  await checkEngine();
  await loadSystems();
  await loadJobs();
}
initializeLocalEngine();
