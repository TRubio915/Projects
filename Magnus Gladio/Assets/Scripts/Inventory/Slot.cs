using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Slot : MonoBehaviour
{

    public GameManager GM;
    public Item slotItem;
    public Image slotImage;
    // public Text slotAmount;
    //public int slotNumber;

    public void Awake()
    {
        GM = GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>();
        //slotAmount = this.GetComponentInChildren<Text>();
    }
    public void ItemOnclicked()
    {

        if (slotItem.itemName.Equals("HealthPotion"))
        {
            GM.IncreasePlayerhealth(20);
        }
        InventoryManager.instance.myBag.itemList.Remove(slotItem);
        InventoryManager.UpdateItemInfo(slotItem.itemInfo);
        Destroy(gameObject);
        InventoryManager.RefreshItem();

    }
}
