// app/api/_diag/ps/route.ts
import { NextResponse } from "next/server";
import { connect } from "@planetscale/database";

export const runtime = "edge";

export async function GET() {
  try {
    const host = process.env.PS_HOST;
    const username = process.env.PS_USERNAME;
    const hasPass = Boolean(process.env.PS_PASSWORD && process.env.PS_PASSWORD!.length >= 4);

    if (!host || !username || !hasPass) {
      return NextResponse.json(
        {
          ok: false,
          step: "env",
          // don't leak secrets; just show presence
          hasHost: Boolean(host),
          hasUsername: Boolean(username),
          hasPassword: hasPass,
        },
        { status: 500 }
      );
    }

    const conn = connect({
      host,
      username,
      password: process.env.PS_PASSWORD!,
    });

    const r = await conn.execute("SELECT 1 AS one");
    return NextResponse.json({ ok: true, result: r.rows });
  } catch (e) {
    const msg =
      e && typeof e === "object" && "message" in e
        ? String((e as { message?: unknown }).message)
        : String(e);
    return NextResponse.json({ ok: false, step: "query", error: msg }, { status: 500 });
  }
}
