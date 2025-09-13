// app/api/tree-image/route.ts
import { NextRequest, NextResponse } from "next/server";
import crypto from "crypto";
import { query } from "@/lib/db";

export const runtime = "nodejs";

/** Copy bytes into a brand-new ArrayBuffer (never SharedArrayBuffer). */
function toSafeArrayBuffer(x: Buffer | Uint8Array): ArrayBuffer {
  const src = Buffer.isBuffer(x)
    ? new Uint8Array(x.buffer, x.byteOffset, x.byteLength)
    : new Uint8Array(x.buffer, x.byteOffset, x.byteLength);
  const ab = new ArrayBuffer(src.byteLength);
  new Uint8Array(ab).set(src);
  return ab;
}

/** Normalize DB binary -> Blob (acceptable BodyInit). */
function toBlob(x: Buffer | Uint8Array, mime: string): Blob {
  const ab = toSafeArrayBuffer(x);
  return new Blob([ab], { type: mime || "application/octet-stream" });
}

// ---------- GET /api/tree-image?job=...&label=final OR ?current=1 ----------
export async function GET(req: NextRequest) {
  try {
    const { searchParams } = new URL(req.url);
    const job = (searchParams.get("job") ?? searchParams.get("job_id") ?? "").trim();
    const label = (searchParams.get("label") ?? "").trim();
    const current = (searchParams.get("current") ?? "").trim() === "1";

    if (!job) {
      return NextResponse.json({ error: "Missing job" }, { status: 400 });
    }

    let rows: Array<{
      image_png: Buffer | Uint8Array | null;
      image_mime: string | null;
    }> = [];

    if (current) {
      const t = await query<{ id: number; label: string | null }>(
        `SELECT id, label
           FROM emtr_trees
          WHERE job_id = ? AND is_current = 1
          ORDER BY id DESC
          LIMIT 1`,
        [job]
      );

      if (t.length) {
        const treeId = t[0].id;
        const lbl = t[0].label || "final";
        rows = await query<{ image_png: Buffer | Uint8Array | null; image_mime: string | null }>(
          `SELECT image_png, image_mime
             FROM emtr_tree_images
            WHERE (tree_id = ?)
               OR (job_id = ? AND label = ?)
            ORDER BY id DESC
            LIMIT 1`,
          [treeId, job, lbl]
        );
      }
    } else if (label) {
      rows = await query<{ image_png: Buffer | Uint8Array | null; image_mime: string | null }>(
        `SELECT image_png, image_mime
           FROM emtr_tree_images
          WHERE job_id = ? AND label = ?
          ORDER BY id DESC
          LIMIT 1`,
        [job, label]
      );
    } else {
      rows = await query<{ image_png: Buffer | Uint8Array | null; image_mime: string | null }>(
        `SELECT image_png, image_mime
           FROM emtr_tree_images
          WHERE job_id = ?
          ORDER BY id DESC
          LIMIT 1`,
        [job]
      );
    }

    if (!rows.length || !rows[0].image_png) {
      return NextResponse.json({ error: "Image not found" }, { status: 404 });
    }

    const mime = rows[0].image_mime || "image/png";
    const body = toBlob(rows[0].image_png, mime);

    return new NextResponse(body, {
      status: 200,
      headers: {
        "content-type": mime,
        "cache-control": "no-store",
      },
    });
  } catch (e) {
    const msg = e instanceof Error ? e.message : String(e);
    return NextResponse.json({ error: msg }, { status: 500 });
  }
}

// ---------- POST /api/tree-image  (multipart/form-data) ----------
// fields: job_id (string, required), label (string, optional; default "final"), image (file, required)
export async function POST(req: NextRequest) {
  try {
    const ct = (req.headers.get("content-type") || "").toLowerCase();
    if (!ct.includes("multipart/form-data")) {
      return NextResponse.json({ ok: false, error: "Expected multipart/form-data" }, { status: 400 });
    }

    const form = await req.formData();
    const job_id = String(form.get("job_id") || "").trim();
    const label = String(form.get("label") || "final").trim() || "final";
    const file = form.get("image");

    if (!job_id) {
      return NextResponse.json({ ok: false, error: "Missing job_id" }, { status: 400 });
    }
    if (!(file instanceof File)) {
      return NextResponse.json({ ok: false, error: "Missing file field 'image'" }, { status: 400 });
    }

    const mime = file.type || "image/png";

    // Image bytes -> base64 string (SQL-safe)
    const imgArrayBuffer = await file.arrayBuffer();
    const imgBase64 = Buffer.from(imgArrayBuffer).toString("base64");

    // SHA-256 -> hex string (SQL-safe)
    const shaHex = crypto.createHash("sha256").update(Buffer.from(imgArrayBuffer)).digest("hex");

    // Optional: link to the current tree row
    const t = await query<{ id: number }>(
      `SELECT id FROM emtr_trees WHERE job_id = ? AND is_current = 1 ORDER BY id DESC LIMIT 1`,
      [job_id]
    );
    const treeId: number | null = t.length ? t[0].id : null;

    // Store as strings; DB reconstructs binary via UNHEX / FROM_BASE64
    await query(
      `INSERT INTO emtr_tree_images (job_id, label, tree_id, image_mime, image_sha256, image_png)
       VALUES (?, ?, ?, ?, UNHEX(?), FROM_BASE64(?))`,
      [job_id, label, treeId, mime, shaHex, imgBase64]
    );

    return NextResponse.json({ ok: true, job_id, label, bytes: Buffer.byteLength(imgBase64, "base64") }, { status: 200 });
  } catch (e) {
    const msg = e instanceof Error ? e.message : String(e);
    return NextResponse.json({ ok: false, error: msg }, { status: 400 });
  }
}
