using System;
using System.IO;
using Microsoft.VisualBasic.FileIO;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Globalization;
using System.Threading.Tasks;

namespace sgRNA
{
	class Program
	{
		static void Main(string[] args)
		{
			Console.WriteLine("======================================================");
			Console.WriteLine("  CRISPR/Cas9用ツール");
			Console.WriteLine("    Copyright (c) 2018–2021 Genta Ito");
			Console.WriteLine("    Version 3.0");
			Console.WriteLine("======================================================");
			Console.WriteLine("  入力ファイル書式：");
			Console.WriteLine("    「標的配列名称+タブ+標的配列（5' -> 3'）+改行」");
			Console.WriteLine("    リストの上下や行間に余計な空白がないように");
			Console.WriteLine("======================================================");
			Console.WriteLine("  出力ファイル：");
			Console.WriteLine("  ・発注用オリゴDNAリストtxtファイル");
			Console.WriteLine("  ・それぞれのオリゴDNA配列SnapGeneファイル");
			Console.WriteLine("  ・アニーリングしたオリゴDNA配列SnapGeneファイル");
			Console.WriteLine("  ・LigationしたベクターマップSnapGeneファイル");
			Console.WriteLine("");
			Console.WriteLine("  入力ファイルがあるフォルダに入力ファイルと同じ名前の");
			Console.WriteLine("  サブフォルダが作成されます");
			Console.WriteLine("======================================================");
			Console.WriteLine("");

			// ドラッグアンドドロップされたファイルのファイルパスを取得
			// 先頭に格納される実行ファイル名を除く
			string[] filePath = Environment.GetCommandLineArgs();
			int startIndex = 0;
			int numberOfFiles = 0;
			for (int i = 0; i < filePath.Length; i++)
			{
				int len = filePath[i].Length;
				if (filePath[i].Substring(len - 3, 3) != "exe")
				{
					startIndex = i;
					break;
				}
			}
			numberOfFiles = filePath.Length - startIndex;

			// エラー処理
			if (numberOfFiles != 1)
			{
				Console.WriteLine("  エラー：ドロップするファイルは1個だけにしてください。");
				Console.WriteLine("  終了するには何かキーを押してください。");
				Console.ReadKey();
				Environment.Exit(0);
			}

			// クエリディクショナリ
			var listNameSeq = new Dictionary<string, string>();
			using (TextFieldParser parser = new TextFieldParser(filePath[startIndex], Encoding.GetEncoding("Shift_JIS")))
			{
				try
				{
					parser.TextFieldType = FieldType.Delimited;
					parser.SetDelimiters("\t");

					while (parser.EndOfData == false)
					{
						string[] data = parser.ReadFields();
						for (int i = 0; i < (data.Length / 2); i++)
						{
							string key, value;
							key = data[i * 2];
							value = data[i * 2 + 1];
							listNameSeq.Add(key, value);
						}
					}
				}
				catch (Exception ex)
				{
					Console.WriteLine("  " + ex.Message);
					Console.WriteLine("  終了するには何かキーを押してください。");
					Console.ReadKey();
					Environment.Exit(0);
				}
			}

			// 結果ディクショナリ
			Dictionary<string, string> results = new Dictionary<string, string>();
			try
			{
				foreach (KeyValuePair<string, string> pair in listNameSeq)
				{
					Dictionary<string, string> result = new Dictionary<string, string>();
					Console.WriteLine("  処理中 {0}: {1}", pair.Key, pair.Value);
					result = MakeOligoSeq(pair.Key, pair.Value);
					results = results.Concat(result).ToDictionary(x => x.Key, x => x.Value);
					Console.WriteLine("  完了");
				}
			}
			catch (Exception ex)
			{
				Console.WriteLine("  " + ex.Message);
				Console.WriteLine("  終了するには何かキーを押してください。");
				Console.ReadKey();
				Environment.Exit(0);
			}

			using (var sw = new StreamWriter(filePath[startIndex], true))
			{
				CultureInfo ci = CultureInfo.CurrentCulture;
				TextInfo ti = ci.TextInfo;

				sw.WriteLine("");
				sw.WriteLine("");
				sw.WriteLine("HTSエクセル入力用 ==============================");
				foreach (KeyValuePair<string, string> kv in results)
				{
					// エクセル入力用の塩基配列は大文字に変換
					sw.WriteLine("{0}" + "\t" + "{1}", kv.Key, ti.ToUpper(kv.Value));
				}
				sw.WriteLine("SnapGene入力用 (FASTA形式)========================");
				foreach (KeyValuePair<string, string> kv in results)
				{
					sw.WriteLine(">" + kv.Key);
					sw.WriteLine(kv.Value);
				}
			}

			Console.WriteLine("  全ての処理が完了しました。");
			Console.WriteLine("  処理結果は入力ファイルに追記されています。");
			Console.WriteLine("  終了するには何かキーを押してください。");
			Console.ReadKey();
			Environment.Exit(0);
		}

		/// <summary>
		/// 標的配列からオリゴDNAセンス鎖とアンチセンス鎖を設計し、返します
		/// </summary>
		/// <param name="name">標的配列名称</param>
		/// <param name="seq">標的配列（20 bp）</param>
		/// <returns></returns>
		static Dictionary<string, string> MakeOligoSeq(string name, string seq)
		{
			Dictionary<string, string> res = new Dictionary<string, string>();

			string name_s, name_as;
			name_s = name + "_s";
			name_as = name + "_as";

			// 標的配列は20塩基でなければならない
			char[] seqArray = seq.ToCharArray();
			if (seqArray.Length != 20)
			{
				throw new ArgumentOutOfRangeException();
			}

			// 標的配列の逆相補鎖の作成（全て大文字）
			char[] seqRevCompArray = new char[seqArray.Length];
			char[] seqArrayCap = new char[seqArray.Length];
			for (var i = 0; i < seqArray.Length; i++)
			{
				if ((seqArray[i] == 'G') || (seqArray[i] == 'g'))
				{
					seqRevCompArray[seqArray.Length - i - 1] = 'C';
					seqArrayCap[i] = 'G';
				}
				else if ((seqArray[i] == 'C') || (seqArray[i] == 'c'))
				{
					seqRevCompArray[seqArray.Length - i - 1] = 'G';
					seqArrayCap[i] = 'C';
				}
				else if ((seqArray[i] == 'A') || (seqArray[i] == 'a'))
				{
					seqRevCompArray[seqArray.Length - i - 1] = 'T';
					seqArrayCap[i] = 'A';
				}
				else if ((seqArray[i] == 'T') || (seqArray[i] == 't'))
				{
					seqRevCompArray[seqArray.Length - i - 1] = 'A';
					seqArrayCap[i] = 'T';
				}
				else
				{
					throw new ArgumentOutOfRangeException();
				}
			}

			// センス鎖の設計（追加されるcaccオーバーハングとgは小文字）
			string seqRes, seqRevCompRes;
			if (seqArrayCap[0] != 'G')
			{
				char[] seqArrayG = new char[21];
				seqArrayG[0] = 'g';
				seqArrayCap.CopyTo(seqArrayG, 1);
				seqRes = new string(seqArrayG);
			}
			else
			{
				seqRes = new string(seqArray);
			}
			res.Add(name_s, "cacc" + seqRes);

			// アンチセンス鎖の設計（追加されるaaacオーバーハングとcは小文字）
			if (seqRevCompArray[19] != 'C')
			{
				char[] seqRevCompArrayG = new char[21];
				seqRevCompArray.CopyTo(seqRevCompArrayG, 0);
				seqRevCompArrayG[20] = 'c';
				seqRevCompRes = new string(seqRevCompArrayG);
			}
			else
			{
				seqRevCompRes = new string(seqRevCompArray);
			}
			res.Add(name_as, "aaac" + seqRevCompRes);

			return res;
		}
	}
}
