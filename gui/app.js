const state = { analyses: [], selected: null, jobs: [], poller: null };
const $ = selector => document.querySelector(selector);
const $$ = selector => [...document.querySelectorAll(selector)];

async function api(path, options = {}) {
  const response = await fetch(path, { headers: { "Content-Type": "application/json" }, ...options });
  const data = await response.json();
  if (!response.ok) throw new Error(data.error || `Request failed (${response.status})`);
  return data;
}

function escapeHtml(value) {
  return String(value ?? "").replace(/[&<>"']/g, c => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" })[c]);
}

let toastTimer;
function toast(title, detail) {
  const node = $("#toast"); node.querySelector("strong").textContent = title; node.querySelector("span").textContent = detail;
  node.classList.add("show"); clearTimeout(toastTimer); toastTimer = setTimeout(() => node.classList.remove("show"), 3600);
}

function switchView(name) {
  $$(".view").forEach(node => node.classList.toggle("active", node.id === `${name}-view`));
  $$(".nav-button").forEach(node => node.classList.toggle("active", node.dataset.view === name));
  const labels = { workbench: "Analysis workbench", jobs: "Execution history", results: "Generated results", setup: "Setup and portability" };
  $("#viewTitle").textContent = labels[name]; document.querySelector(".sidebar").classList.remove("open"); $("#mobileShade").classList.remove("open");
  if (name === "jobs" || name === "results") loadJobs();
}

function renderStatus(status) {
  $("#repositoryPath").textContent = status.repository;
  const engine = $("#engineStatus"); engine.classList.add("online"); engine.querySelector("strong").textContent = "Engine ready";
  engine.querySelector("small").textContent = `Python ${status.dependencies.python}`;
  const labels = [["ml_stack", "ML stack"], ["vmd", "VMD"], ["namd2", "NAMD"]];
  $("#dependencyStrip").innerHTML = labels.map(([key, label]) => `<span class="${status.dependencies[key] ? "ready" : "missing"}">${label}: ${status.dependencies[key] ? "ready" : "not found"}</span>`).join("");
}

function filteredAnalyses() {
  const query = $("#analysisSearch").value.toLowerCase(); const group = $("#groupFilter").value;
  return state.analyses.filter(item => (!group || item.group === group) && (!query || `${item.name} ${item.summary} ${item.group}`.toLowerCase().includes(query)));
}

function renderCatalog() {
  const items = filteredAnalyses();
  $("#analysisList").innerHTML = items.length ? items.map((item, index) => `<button class="analysis-card ${state.selected?.id === item.id ? "active" : ""}" data-analysis="${item.id}"><span class="analysis-index">${String(index + 1).padStart(2, "0")}</span><span><strong>${escapeHtml(item.name)}</strong><small>${escapeHtml(item.summary)}</small></span><em>${escapeHtml(item.group)}</em></button>`).join("") : `<div class="empty-state"><strong>No matching workflows</strong><span>Adjust the search or group filter.</span></div>`;
  $$(".analysis-card").forEach(button => button.addEventListener("click", () => selectAnalysis(button.dataset.analysis)));
}

function inputMarkup(spec) {
  const pathClass = ["directory", "file"].includes(spec.type) ? "path-field" : "";
  const kind = ["integer", "number"].includes(spec.type) ? "number" : "text";
  const attributes = [`data-key="${spec.key}"`, `type="${kind}"`, `value="${escapeHtml(spec.default)}"`];
  if (spec.min !== undefined) attributes.push(`min="${spec.min}"`); if (spec.max !== undefined) attributes.push(`max="${spec.max}"`);
  if (spec.type === "number") attributes.push("step=\"any\"");
  return `<label class="field ${pathClass}"><span>${escapeHtml(spec.label)}<i>${escapeHtml(spec.unit || spec.type)}</i></span><input ${attributes.join(" ")} ${spec.optional ? "" : "required"}></label>`;
}

function selectAnalysis(id) {
  state.selected = state.analyses.find(item => item.id === id);
  $("#selectedGroup").textContent = state.selected.group; $("#selectedName").textContent = state.selected.name;
  $("#selectedSummary").textContent = state.selected.summary; $("#expectedOutputs").textContent = state.selected.outputs;
  $("#formFields").innerHTML = state.selected.fields.map(inputMarkup).join(""); $("#runButton").disabled = false;
  $$("#formFields input").forEach(input => input.addEventListener("input", updatePreview));
  renderCatalog(); updatePreview();
}

function formValues() {
  return Object.fromEntries($$("#formFields input").map(input => [input.dataset.key, input.value]));
}

function updatePreview() {
  if (!state.selected) return;
  const parts = [`python3 src/${state.selected.script}`];
  for (const [key, value] of Object.entries(formValues())) if (value !== "") parts.push(`--${key} ${/\s/.test(value) ? `"${value}"` : value}`);
  parts.push("--output-dir gui/local-results/<job-id>"); $("#commandPreview").textContent = parts.join(" ");
}

function validateClient() {
  if (!state.selected) throw new Error("Select an analysis first.");
  let valid = true;
  $$("#formFields input").forEach(input => { const bad = input.required && !input.value.trim(); input.classList.toggle("invalid", bad); valid = valid && !bad; });
  if (!valid) throw new Error("Complete the required fields.");
  return formValues();
}

async function runAnalysis() {
  try {
    const values = validateClient(); const button = $("#runButton"); button.disabled = true; button.textContent = "Starting...";
    const job = await api("/api/run", { method: "POST", body: JSON.stringify({ analysis: state.selected.id, values }) });
    toast("Analysis queued", `${job.name} is running in ${job.id}.`); await loadJobs(); switchView("jobs");
  } catch (error) { toast("Unable to start", error.message); }
  finally { $("#runButton").disabled = false; $("#runButton").textContent = "Run analysis"; }
}

function formatDate(value) { return value ? new Date(value).toLocaleString([], { dateStyle: "medium", timeStyle: "short" }) : "—"; }
function formatBytes(bytes) { if (bytes < 1024) return `${bytes} B`; if (bytes < 1048576) return `${(bytes / 1024).toFixed(1)} KB`; return `${(bytes / 1048576).toFixed(1)} MB`; }

function renderJobs() {
  const active = state.jobs.filter(job => ["queued", "running"].includes(job.status)).length; $("#activeJobCount").textContent = active;
  $("#jobList").innerHTML = state.jobs.length ? state.jobs.map(job => `<div class="job-row"><div><strong>${escapeHtml(job.name)}</strong><small>${escapeHtml(job.id)}</small></div><span>${formatDate(job.created_at)}</span><span class="status ${job.status}"><i></i>${escapeHtml(job.status)}</span><span>${escapeHtml(job.message)}</span><button data-log="${job.id}">Inspect</button></div>`).join("") : `<div class="empty-state"><strong>No jobs yet</strong><span>Configure a workflow in the workbench.</span></div>`;
  $$('[data-log]').forEach(button => button.addEventListener("click", () => showLog(button.dataset.log)));
  renderResults();
  if (active && !state.poller) state.poller = setInterval(loadJobs, 1800);
  if (!active && state.poller) { clearInterval(state.poller); state.poller = null; }
}

async function loadJobs() { try { state.jobs = (await api("/api/jobs")).jobs || []; renderJobs(); } catch (error) { toast("Job service unavailable", error.message); } }

function showLog(id) {
  const job = state.jobs.find(item => item.id === id); if (!job) return;
  $("#logPanel").hidden = false; $("#logTitle").textContent = `${job.name} · ${job.status}`;
  $("#logContent").textContent = `$ ${job.command.map(item => /\s/.test(item) ? `"${item}"` : item).join(" ")}\n\n${job.log_tail || "Waiting for log output..."}`;
  $("#logPanel").scrollIntoView({ behavior: "smooth", block: "nearest" });
}

function renderResults() {
  const complete = state.jobs.filter(job => job.status === "complete" && job.outputs?.length);
  $("#resultBrowser").innerHTML = complete.length ? complete.map(job => `<section class="result-group"><header><h2>${escapeHtml(job.name)}</h2><span>${job.outputs.length} files · ${formatDate(job.finished_at)}</span></header>${job.outputs.map(file => `<div class="file-row"><code title="${escapeHtml(file.path)}">${escapeHtml(file.path)}</code><span>${formatBytes(file.bytes)}</span><a class="file-link" href="/download/${job.id}/${encodeURI(file.path)}">Download</a></div>`).join("")}</section>`).join("") : `<div class="empty-state"><strong>No completed jobs yet</strong><span>Run an analysis to populate this view.</span></div>`;
}

function resetDefaults() { if (state.selected) selectAnalysis(state.selected.id); toast("Defaults restored", "Package paths and recommended parameters were restored."); }

async function init() {
  try {
    const [status, catalog] = await Promise.all([api("/api/status"), api("/api/analyses")]); renderStatus(status); state.analyses = catalog.analyses;
    const groups = [...new Set(state.analyses.map(item => item.group))]; $("#groupFilter").innerHTML += groups.map(group => `<option>${escapeHtml(group)}</option>`).join("");
    renderCatalog(); selectAnalysis(state.analyses[0].id); await loadJobs();
  } catch (error) { $("#engineStatus strong").textContent = "Engine offline"; $("#engineStatus small").textContent = "Start gui/server.py"; toast("Local engine unavailable", error.message); }
}

$$('.nav-button').forEach(button => button.addEventListener("click", () => switchView(button.dataset.view)));
$("#analysisSearch").addEventListener("input", renderCatalog); $("#groupFilter").addEventListener("change", renderCatalog);
$("#runButton").addEventListener("click", runAnalysis); $("#validateButton").addEventListener("click", () => { try { validateClient(); toast("Inputs look complete", "The server will verify local paths when the job starts."); } catch (error) { toast("Review inputs", error.message); } });
$("#useDefaultsButton").addEventListener("click", resetDefaults); $("#refreshJobsButton").addEventListener("click", loadJobs); $("#closeLogButton").addEventListener("click", () => $("#logPanel").hidden = true);
$("#menuButton").addEventListener("click", () => { document.querySelector(".sidebar").classList.add("open"); $("#mobileShade").classList.add("open"); });
$("#mobileShade").addEventListener("click", () => { document.querySelector(".sidebar").classList.remove("open"); $("#mobileShade").classList.remove("open"); });
init();
