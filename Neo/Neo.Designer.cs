namespace Neo
{
    partial class Neo
    {
        /// <summary>
        /// Required designer variable.
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Clean up any resources being used.
        /// </summary>
        /// <param name="disposing">true if managed resources should be disposed; otherwise, false.</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form Designer generated code

        /// <summary>
        /// Required method for Designer support - do not modify
        /// the contents of this method with the code editor.
        /// </summary>
        private void InitializeComponent()
        {
            this.bButton = new System.Windows.Forms.Button();
            this.FASTAButton = new System.Windows.Forms.Button();
            this.yButton = new System.Windows.Forms.Button();
            this.btxtBox = new System.Windows.Forms.TextBox();
            this.ytxtBox = new System.Windows.Forms.TextBox();
            this.FASTAtxtBox = new System.Windows.Forms.TextBox();
            this.label2 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.StartButton = new System.Windows.Forms.Button();
            this.label1 = new System.Windows.Forms.Label();
            this.progressBar1 = new System.Windows.Forms.ProgressBar();
            this.StatustxtBox = new System.Windows.Forms.TextBox();
            this.label4 = new System.Windows.Forms.Label();
            this.backgroundWorker1 = new System.ComponentModel.BackgroundWorker();
            this.MS1Tolerance = new System.Windows.Forms.TextBox();
            this.MS2Tolerance = new System.Windows.Forms.TextBox();
            this.label5 = new System.Windows.Forms.Label();
            this.label6 = new System.Windows.Forms.Label();
            this.label11 = new System.Windows.Forms.Label();
            this.bCheckBox = new System.Windows.Forms.CheckBox();
            this.cCheckBox = new System.Windows.Forms.CheckBox();
            this.yCheckBox = new System.Windows.Forms.CheckBox();
            this.zdotCheckBox = new System.Windows.Forms.CheckBox();
            this.checkBox1 = new System.Windows.Forms.CheckBox();
            this.SuspendLayout();
            // 
            // bButton
            // 
            this.bButton.Location = new System.Drawing.Point(465, 82);
            this.bButton.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.bButton.Name = "bButton";
            this.bButton.Size = new System.Drawing.Size(51, 20);
            this.bButton.TabIndex = 0;
            this.bButton.Text = "Browse";
            this.bButton.UseVisualStyleBackColor = true;
            this.bButton.Click += new System.EventHandler(this.bButton_Click);
            // 
            // FASTAButton
            // 
            this.FASTAButton.Location = new System.Drawing.Point(465, 206);
            this.FASTAButton.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.FASTAButton.Name = "FASTAButton";
            this.FASTAButton.Size = new System.Drawing.Size(51, 19);
            this.FASTAButton.TabIndex = 1;
            this.FASTAButton.Text = "Browse";
            this.FASTAButton.UseVisualStyleBackColor = true;
            this.FASTAButton.Click += new System.EventHandler(this.FASTAButton_Click);
            // 
            // yButton
            // 
            this.yButton.Location = new System.Drawing.Point(465, 129);
            this.yButton.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.yButton.Name = "yButton";
            this.yButton.Size = new System.Drawing.Size(51, 19);
            this.yButton.TabIndex = 2;
            this.yButton.Text = "Browse";
            this.yButton.UseVisualStyleBackColor = true;
            this.yButton.Click += new System.EventHandler(this.yButton_Click);
            // 
            // btxtBox
            // 
            this.btxtBox.Location = new System.Drawing.Point(51, 82);
            this.btxtBox.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.btxtBox.Name = "btxtBox";
            this.btxtBox.Size = new System.Drawing.Size(387, 20);
            this.btxtBox.TabIndex = 3;
            this.btxtBox.TextChanged += new System.EventHandler(this.btxtBox_TextChanged);
            // 
            // ytxtBox
            // 
            this.ytxtBox.Location = new System.Drawing.Point(52, 130);
            this.ytxtBox.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.ytxtBox.Name = "ytxtBox";
            this.ytxtBox.Size = new System.Drawing.Size(386, 20);
            this.ytxtBox.TabIndex = 4;
            this.ytxtBox.TextChanged += new System.EventHandler(this.ytxtBox_TextChanged);
            // 
            // FASTAtxtBox
            // 
            this.FASTAtxtBox.Location = new System.Drawing.Point(51, 207);
            this.FASTAtxtBox.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.FASTAtxtBox.Name = "FASTAtxtBox";
            this.FASTAtxtBox.Size = new System.Drawing.Size(387, 20);
            this.FASTAtxtBox.TabIndex = 5;
            this.FASTAtxtBox.TextChanged += new System.EventHandler(this.FASTAtxtBox_TextChanged);
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(50, 114);
            this.label2.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(76, 13);
            this.label2.TabIndex = 7;
            this.label2.Text = "C Terminal File";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(50, 191);
            this.label3.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(60, 13);
            this.label3.TabIndex = 8;
            this.label3.Text = "FASTA File";
            // 
            // StartButton
            // 
            this.StartButton.Location = new System.Drawing.Point(265, 344);
            this.StartButton.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.StartButton.Name = "StartButton";
            this.StartButton.Size = new System.Drawing.Size(90, 37);
            this.StartButton.TabIndex = 9;
            this.StartButton.Text = "Analyze";
            this.StartButton.UseVisualStyleBackColor = true;
            this.StartButton.Click += new System.EventHandler(this.StartButton_Click);
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Location = new System.Drawing.Point(50, 66);
            this.label1.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(77, 13);
            this.label1.TabIndex = 10;
            this.label1.Text = "N Terminal File";
            this.label1.Click += new System.EventHandler(this.label1_Click);
            // 
            // progressBar1
            // 
            this.progressBar1.Location = new System.Drawing.Point(393, 416);
            this.progressBar1.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.progressBar1.Name = "progressBar1";
            this.progressBar1.Size = new System.Drawing.Size(123, 18);
            this.progressBar1.TabIndex = 11;
            this.progressBar1.Click += new System.EventHandler(this.progressBar1_Click);
            // 
            // StatustxtBox
            // 
            this.StatustxtBox.Location = new System.Drawing.Point(51, 416);
            this.StatustxtBox.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.StatustxtBox.Name = "StatustxtBox";
            this.StatustxtBox.Size = new System.Drawing.Size(305, 20);
            this.StatustxtBox.TabIndex = 12;
            this.StatustxtBox.TextChanged += new System.EventHandler(this.StatustxtBox_TextChanged);
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(46, 400);
            this.label4.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(63, 13);
            this.label4.TabIndex = 13;
            this.label4.Text = "Run Status:";
            this.label4.Click += new System.EventHandler(this.label4_Click);
            // 
            // backgroundWorker1
            // 
            this.backgroundWorker1.DoWork += new System.ComponentModel.DoWorkEventHandler(this.backgroundWorker1_DoWork);
            // 
            // MS1Tolerance
            // 
            this.MS1Tolerance.Location = new System.Drawing.Point(188, 239);
            this.MS1Tolerance.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.MS1Tolerance.Name = "MS1Tolerance";
            this.MS1Tolerance.Size = new System.Drawing.Size(43, 20);
            this.MS1Tolerance.TabIndex = 14;
            this.MS1Tolerance.TextChanged += new System.EventHandler(this.MS1Tolerance_TextChanged);
            // 
            // MS2Tolerance
            // 
            this.MS2Tolerance.Location = new System.Drawing.Point(188, 266);
            this.MS2Tolerance.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.MS2Tolerance.Name = "MS2Tolerance";
            this.MS2Tolerance.Size = new System.Drawing.Size(43, 20);
            this.MS2Tolerance.TabIndex = 15;
            this.MS2Tolerance.TextChanged += new System.EventHandler(this.MS2Tolerance_TextChanged);
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(49, 240);
            this.label5.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(135, 13);
            this.label5.TabIndex = 16;
            this.label5.Text = "Precursor Tolerance (ppm):";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(50, 267);
            this.label6.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(128, 13);
            this.label6.TabIndex = 17;
            this.label6.Text = "Fragment Tolerance (Da):";
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(278, 240);
            this.label11.Margin = new System.Windows.Forms.Padding(2, 0, 2, 0);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(79, 13);
            this.label11.TabIndex = 22;
            this.label11.Text = "Ions Searched:";
            // 
            // bCheckBox
            // 
            this.bCheckBox.AutoSize = true;
            this.bCheckBox.Checked = true;
            this.bCheckBox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.bCheckBox.Location = new System.Drawing.Point(362, 240);
            this.bCheckBox.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.bCheckBox.Name = "bCheckBox";
            this.bCheckBox.Size = new System.Drawing.Size(32, 17);
            this.bCheckBox.TabIndex = 23;
            this.bCheckBox.Text = "b";
            this.bCheckBox.UseVisualStyleBackColor = true;
            this.bCheckBox.CheckedChanged += new System.EventHandler(this.bCheckBox_CheckedChanged);
            // 
            // cCheckBox
            // 
            this.cCheckBox.AutoSize = true;
            this.cCheckBox.Location = new System.Drawing.Point(362, 264);
            this.cCheckBox.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.cCheckBox.Name = "cCheckBox";
            this.cCheckBox.Size = new System.Drawing.Size(32, 17);
            this.cCheckBox.TabIndex = 24;
            this.cCheckBox.Text = "c";
            this.cCheckBox.UseVisualStyleBackColor = true;
            this.cCheckBox.CheckedChanged += new System.EventHandler(this.cCheckBox_CheckedChanged);
            // 
            // yCheckBox
            // 
            this.yCheckBox.AutoSize = true;
            this.yCheckBox.Checked = true;
            this.yCheckBox.CheckState = System.Windows.Forms.CheckState.Checked;
            this.yCheckBox.Location = new System.Drawing.Point(394, 240);
            this.yCheckBox.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.yCheckBox.Name = "yCheckBox";
            this.yCheckBox.Size = new System.Drawing.Size(31, 17);
            this.yCheckBox.TabIndex = 25;
            this.yCheckBox.Text = "y";
            this.yCheckBox.UseVisualStyleBackColor = true;
            this.yCheckBox.CheckedChanged += new System.EventHandler(this.yCheckBox_CheckedChanged);
            // 
            // zdotCheckBox
            // 
            this.zdotCheckBox.AutoSize = true;
            this.zdotCheckBox.Location = new System.Drawing.Point(394, 264);
            this.zdotCheckBox.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.zdotCheckBox.Name = "zdotCheckBox";
            this.zdotCheckBox.Size = new System.Drawing.Size(46, 17);
            this.zdotCheckBox.TabIndex = 26;
            this.zdotCheckBox.Text = "zdot";
            this.zdotCheckBox.UseVisualStyleBackColor = true;
            this.zdotCheckBox.CheckedChanged += new System.EventHandler(this.zdotCheckBox_CheckedChanged);
            // 
            // checkBox1
            // 
            this.checkBox1.AutoSize = true;
            this.checkBox1.Location = new System.Drawing.Point(53, 291);
            this.checkBox1.Margin = new System.Windows.Forms.Padding(2);
            this.checkBox1.Name = "checkBox1";
            this.checkBox1.Size = new System.Drawing.Size(94, 17);
            this.checkBox1.TabIndex = 27;
            this.checkBox1.Text = "Decoy Search";
            this.checkBox1.UseVisualStyleBackColor = true;
            this.checkBox1.CheckedChanged += new System.EventHandler(this.checkBox1_CheckedChanged);
            // 
            // Neo
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(572, 519);
            this.Controls.Add(this.checkBox1);
            this.Controls.Add(this.zdotCheckBox);
            this.Controls.Add(this.yCheckBox);
            this.Controls.Add(this.cCheckBox);
            this.Controls.Add(this.bCheckBox);
            this.Controls.Add(this.label11);
            this.Controls.Add(this.label6);
            this.Controls.Add(this.label5);
            this.Controls.Add(this.MS2Tolerance);
            this.Controls.Add(this.MS1Tolerance);
            this.Controls.Add(this.label4);
            this.Controls.Add(this.StatustxtBox);
            this.Controls.Add(this.progressBar1);
            this.Controls.Add(this.label1);
            this.Controls.Add(this.StartButton);
            this.Controls.Add(this.label3);
            this.Controls.Add(this.label2);
            this.Controls.Add(this.FASTAtxtBox);
            this.Controls.Add(this.ytxtBox);
            this.Controls.Add(this.btxtBox);
            this.Controls.Add(this.yButton);
            this.Controls.Add(this.FASTAButton);
            this.Controls.Add(this.bButton);
            this.Margin = new System.Windows.Forms.Padding(2, 2, 2, 2);
            this.Name = "Neo";
            this.Text = "Neo - Fusion Peptide Software";
            this.Load += new System.EventHandler(this.Neo_Load);
            this.ResumeLayout(false);
            this.PerformLayout();

        }

        #endregion

        private System.Windows.Forms.Button bButton;
        private System.Windows.Forms.Button FASTAButton;
        private System.Windows.Forms.Button yButton;
        private System.Windows.Forms.TextBox btxtBox;
        private System.Windows.Forms.TextBox ytxtBox;
        private System.Windows.Forms.TextBox FASTAtxtBox;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.Button StartButton;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.ProgressBar progressBar1;
        private System.Windows.Forms.TextBox StatustxtBox;
        private System.Windows.Forms.Label label4;
        public System.ComponentModel.BackgroundWorker backgroundWorker1;
        private System.Windows.Forms.TextBox MS1Tolerance;
        private System.Windows.Forms.TextBox MS2Tolerance;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.CheckBox bCheckBox;
        private System.Windows.Forms.CheckBox cCheckBox;
        private System.Windows.Forms.CheckBox yCheckBox;
        private System.Windows.Forms.CheckBox zdotCheckBox;
        private System.Windows.Forms.CheckBox checkBox1;
    }
}

