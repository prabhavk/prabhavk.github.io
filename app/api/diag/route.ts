// app/api/diag/route.ts
import { NextResponse } from "next/server";
import { connect } from "@planetscale/database";

export const runtime = "nodejs";

export async function GET() {
  try {
    const { PS_HOST, PS_USERNAME, PS_PASSWORD } = process.env;
    if (!PS_HOST || !PS_USERNAME || !PS_PASSWORD) {
      return NextResponse.json(
        { ok: false, step: "env", hasHost: !!PS_HOST, hasUsername: !!PS_USERNAME, hasPassword: !!PS_PASSWORD },
        { status: 500 }
      );
    }
    const conn = connect({ host: PS_HOST, username: PS_USERNAME, password: PS_PASSWORD });
    const r = await conn.execute("SELECT 1 AS one");
    return NextResponse.json({ ok: true, result: r.rows });
  } catch (e) {
    return NextResponse.json({ ok: false, step: "query", error: String(e) }, { status: 500 });
  }
}
